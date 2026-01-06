package main

import (
	"os"
	"os/signal"
	"fmt"
	"unsafe"
	"strconv"
	"net"
	"net/http"
	"strings"
	"regexp"
	"syscall"
	"time"
)

/*
#cgo LDFLAGS: kalloc.o gfa-base.o gfa-io.o gfa-util.o -lz -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa-priv.h"

char *gfa_extract(gfa_t *g, const char *str, int step, const char *ori_gene)
{
	int i, n, *seg, n_seg, len;
	char **list, *out;
	gfa_t *f;
	list = gfa_read_list(str, &n);
	seg = gfa_list2seg(g, n, list, &n_seg);
	for (i = 0; i < n; ++i) free(list[i]);
	free(list);
	if (step > 0) {
		int32_t n_seg0 = n_seg, *seg0 = seg;
		seg = gfa_sub_extend(g, n_seg0, seg0, step, &n_seg);
		free(seg0);
	}
	f = gfa_subview2(g, n_seg, seg, 1);
	gfa_walk_flip(f, ori_gene);
	out = gfa_write(f, GFA_O_NO_SEQ, &len);
	gfa_subview_destroy(f);
	return out;
}
*/
import "C"

/****************
 * BSD getopt() *
 ****************/

var optind int = 1
var getopt_place int = -1

func getopt(args []string, ostr string) (int, string) {
	if getopt_place == -1 { // update scanning pointer
		if optind >= len(args) || args[optind][0] != '-' {
			getopt_place = -1
			return -1, ""
		}
		if optind < len(args) {
			getopt_place = 0
		}
		if getopt_place + 1 < len(args[optind]) {
			getopt_place += 1
			if args[optind][getopt_place] == '-' { // found "--"
				optind += 1
				getopt_place = -1
				return -1, ""
			}
		}
	}
	optopt := args[optind][getopt_place];
	getopt_place += 1
	oli, arg := strings.IndexByte(ostr, optopt), "";
	if optopt == ':' || oli < 0 {
		if optopt == '-' {
			return -1, ""
		}
		if getopt_place < 0 {
			optind += 1
		}
		return '?', ""
	}
	if oli + 1 >= len(ostr) || ostr[oli+1] != ':' {
		if getopt_place < 0 || getopt_place >= len(args[optind]) {
			optind += 1
			getopt_place = -1
		}
	} else {
		if getopt_place >= 0 && getopt_place < len(args[optind]) {
			arg = args[optind][getopt_place:]
		} else if optind += 1; len(args) <= optind { // no arguments
			getopt_place = -1
			if len(ostr) > 0 && ostr[0] == ':' {
				return ':', ""
			}
			return '?', ""
		} else {
			arg = args[optind]
		}
		getopt_place = -1
		optind += 1
	}
	return int(optopt), arg
}

/************
 * Handlers *
 ************/

var gfa_server_port string = "8000";
var gfa_server_socket string = "";
var gfa_server_ver string = "r294";
var gfa_endpoint string = "/view";
var gfa_graphs map[string]*C.gfa_t;
var gfa_graph_list []string;
var gfa_graph_default *C.gfa_t;
var gfa_js_dir string = "js/";
var gfa_html_dir string = "";
var gfa_max_out int = 200000;

func gfa_print_page(w http.ResponseWriter, r *http.Request, graph_str string) {
	graph, genes, step, ori := "", "", "3", "";
	if len(r.Form["graph"]) > 0 {
		graph = r.Form["graph"][0];
	}
	if len(r.Form["gene"]) > 0 {
		genes = r.Form["gene"][0];
	}
	if len(r.Form["step"]) > 0 {
		step = r.Form["step"][0];
	}
	if len(r.Form["ori"]) > 0 {
		ori = r.Form["ori"][0];
	}
	w.Header().Set("Content-Type", "text/html; charset=utf-8");
	fmt.Fprintln(w, `<title>GFA view</title>`);
	fmt.Fprintln(w, `<style type="text/css">`);
	fmt.Fprintln(w, `  #canvas_graph,#canvas_walk { border: 1px solid #000; }`);
	fmt.Fprintln(w, `  body { font-family: "helvetica neue", helvetica, arial, sans-serif; font-size: 0.9em; }`);
	fmt.Fprintln(w, `  textarea { font-family: monospace; }`);
	fmt.Fprintln(w, `</style>`);
	fmt.Fprintln(w, `<script language="JavaScript" src="js/gfa.js"></script>`);
	fmt.Fprintln(w, `<script language="JavaScript" src="js/gfa-plot.js"></script>`);
	fmt.Fprintln(w, `<body onLoad="plot();">`);
	fmt.Fprintln(w, `<form action="` + gfa_endpoint + `" method="GET">`);
	fmt.Fprintln(w, `  Graph: <select name="graph">`);
	for i := 0; i < len(gfa_graph_list); i++ {
		selected := "";
		if gfa_graph_list[i] == graph {
			selected = " selected";
		}
		fmt.Fprintln(w, `    <option value="` + gfa_graph_list[i] + `"` + selected + `>` + gfa_graph_list[i] + `</option>`);
	}
	fmt.Fprintln(w, `  </select>&nbsp;`);
	fmt.Fprintln(w, `  genes: <input name="gene" size="30" value="` + genes + `"/>&nbsp;`);
	fmt.Fprintln(w, `  neighbors: <input name="step" size="5" value="` + step + `"/>&nbsp;`);
	fmt.Fprintln(w, `  strand: <input name="ori" size="10" value="` + ori + `"/>&nbsp;`);
	fmt.Fprintln(w, `  <input type="submit" value="Retrieve"/>`);
	fmt.Fprintln(w, `</form>`);
	fmt.Fprintln(w, `<hr/>`);
	if graph_str == "" {
		fmt.Fprintln(w, `<h2>Instructions</h2>`);
		fmt.Fprintln(w, `<p>Select a graph, provide one or multiple colocalized genes (<b>EXACT</b> gene names, <b>case-sensitive</b>) and click the`);
		fmt.Fprintln(w, `"Retrieve" button to extract a subgraph around the genes and plot it.`);
		fmt.Fprintln(w, `"Neighbors" controls how many neighboring genes to explore. "Strand"`);
		fmt.Fprintln(w, `orients haplotypes in the given gene strand.</p>`);
		fmt.Fprintln(w, `<p>In the plot, each arrow corresponds to a gene. An empty arrow suggests frameshifts or`);
		fmt.Fprintln(w, `in-frame stop codons, which may be biological or alignment errors. Note that`);
		fmt.Fprintln(w, `inputting genes on different chromosomes or distant apart may lead to`);
		fmt.Fprintln(w, `undesired plots. If you don't see any haplotypes, unclick "filter fragmented contigs".</p>`);
		fmt.Fprintln(w, `<p>Once you see the graph, you may click the "Replot" button to randomize`);
		fmt.Fprintln(w, `node colors. Replotting does not incur server load and is the preferred way`);
		fmt.Fprintln(w, `to adjust plotting.</p>`);
		return;
	}
	fmt.Fprintln(w, `<p>Plot setting: <input type="checkbox" id="merge_walk" checked/>merge identical paths &#xFF5C;`);
	fmt.Fprintln(w, `<input type="checkbox" id="uniq_walk" checked/>filter fragmented contigs &#xFF5C;`);
	fmt.Fprintln(w, `<input type="checkbox" id="rand_color"/>random color &#xFF5C;&nbsp;`);
	fmt.Fprintln(w, `<input type="button" value="Replot" onClick="plot();"></p>`);
	fmt.Fprintln(w, `<p><canvas id="canvas_walk" width="800" height="100"></canvas></p>`);
	fmt.Fprintln(w, `<p><canvas id="canvas_graph" width="800" height="100"></canvas></p>`);
	fmt.Fprintln(w, `<textarea id="gfa-text" readonly rows="15" cols="110">`);
	fmt.Fprintf(w, graph_str);
	fmt.Fprintln(w, `</textarea>`);
	fmt.Fprintln(w, `<script language="JavaScript">`);
	fmt.Fprintln(w, `function plot() {`);
	fmt.Fprintln(w, `  var gfa_conf = gfa_plot_conf();`);
	fmt.Fprintln(w, `  gfa_conf.merge_walk = document.getElementById('merge_walk').checked;`);
	fmt.Fprintln(w, `  gfa_conf.uniq_walk = document.getElementById('uniq_walk').checked;`);
	fmt.Fprintln(w, `  var gfa_text = document.getElementById("gfa-text").value;`);
	fmt.Fprintln(w, `  var gfa_graph = gfa_parse(gfa_text);`);
	fmt.Fprintln(w, `  gfa_assign_color(gfa_graph, document.getElementById('rand_color').checked);`);
	fmt.Fprintln(w, `  gfa_plot_graph(document.getElementById("canvas_graph"), gfa_conf, gfa_graph);`);
	fmt.Fprintln(w, `  gfa_plot_walk(document.getElementById("canvas_walk"), gfa_conf, gfa_graph);`);
	fmt.Fprintln(w, `}`);
	fmt.Fprintln(w, `</script>`);
	fmt.Fprintln(w, `</body>`);
}

func gfa_server_query(w http.ResponseWriter, r *http.Request) {
	r.ParseForm();
	start_time := time.Now().UnixNano();
	fmt.Fprintf(os.Stderr, "[%d] got request: %s\n", start_time, r.Form);
	defer fmt.Fprintf(os.Stderr, "[%d] responded %d\n", time.Now().UnixNano(), start_time);
	if len(r.Form) == 0 {
		gfa_print_page(w, r, "");
		return;
	}
	g := gfa_graph_default;
	step := 3;
	if len(r.Form["graph"]) > 0 { // set graph
		tmp, ok := gfa_graphs[r.Form["graph"][0]];
		if !ok {
			http.Error(w, "400 Bad Request: failed to find graph '" + r.Form["graph"][0] + "'", 400);
			return;
		}
		g = tmp;
	}
	if len(r.Form["step"]) > 0 { // set radius
		i, _ := strconv.Atoi(r.Form["step"][0]);
		if i < 0 {
			http.Error(w, "400 Bad Request: 'r' shouldn't be negative", 400);
			return;
		}
		step = i;
	}
	if len(r.Form["gene"]) > 0 {
		ori := "\t";
		if len(r.Form["ori"]) > 0 {
			ori = r.Form["ori"][0];
		}
		cstr := C.CString(r.Form["gene"][0]);
		cori := C.CString(ori);
		out := C.gfa_extract(g, cstr, C.int(step), cori);
		C.free(unsafe.Pointer(cori));
		C.free(unsafe.Pointer(cstr));
		ret := C.GoString(out);
		C.free(unsafe.Pointer(out));
		if len(ret) > gfa_max_out {
			http.Error(w, "400 Bad Request: subgraph is over 100,000 bytes in size", 400);
			return;
		}
		gfa_print_page(w, r, ret);
	}
}

/*****************
 * Main function *
 *****************/

func main() {
	// set PORT
	if os.Getenv("PORT") != "" {
		gfa_server_port = os.Getenv("PORT");
	}

	// parse command line options
	for {
		opt, arg := getopt(os.Args, "p:s:e:j:d:vm:");
		if opt == 'p' {
			gfa_server_port = arg;
		} else if opt == 's' {
			gfa_server_socket = arg;
		} else if opt == 'e' {
			gfa_endpoint = arg;
		} else if opt == 'j' {
			gfa_js_dir = arg;
		} else if opt == 'd' {
			gfa_html_dir = arg;
		} else if opt == 'v' {
			fmt.Println(gfa_server_ver);
			os.Exit(0);
		} else if opt == 'm' {
			gfa_max_out, _ = strconv.Atoi(arg);
		} else if opt < 0 {
			break;
		}
	}
	if optind == len(os.Args) {
		fmt.Fprintln(os.Stderr, "Usage: gfa-server [options] <graph1.gfa> [graph2.gfa [...]]");
		fmt.Fprintln(os.Stderr, "Options:");
		fmt.Fprintf(os.Stderr, "  -p INT    port number [%s or from $PORT env]\n", gfa_server_port);
		fmt.Fprintf(os.Stderr, "  -s STR    path to Unix domain socket to use instead of a port\n");
		fmt.Fprintf(os.Stderr, "  -j DIR    directory to gfa javascript files [%s]\n", gfa_js_dir);
		fmt.Fprintf(os.Stderr, "  -d DIR    directory to HTML pages to be served at \"/\" []\n");
		fmt.Fprintf(os.Stderr, "  -e STR    endpoint [%s]\n", gfa_endpoint);
		fmt.Fprintf(os.Stderr, "  -m INT    max GFA output size [%d]\n", gfa_max_out);
		fmt.Fprintf(os.Stderr, "  -v        print version number\n");
		os.Exit(1);
	}

	re_base := regexp.MustCompile(`(.*/)?([^/]+)$`);
	re_gz := regexp.MustCompile(`\.gz$`);
	re_gfa := regexp.MustCompile(`\.gfa$`);
	gfa_graphs = make(map[string]*C.gfa_t)
	for i := optind; i < len(os.Args); i++ {
		fn := C.CString(os.Args[i]);
		defer C.free(unsafe.Pointer(fn));
		key := re_base.ReplaceAllString(os.Args[i], "$2");
		key = re_gz.ReplaceAllString(key, "");
		key = re_gfa.ReplaceAllString(key, "");
		fmt.Fprintf(os.Stderr, "[%d] read graph '%s' from '%s'\n", time.Now().UnixNano(), key, os.Args[i]);
		g := C.gfa_read(fn);
		defer C.gfa_destroy(g);
		gfa_graphs[key] = g; // TODO: check duplicate names!!!
		gfa_graph_list = append(gfa_graph_list, key);
		if i == optind {
			gfa_graph_default = g;
		}
	}

	http.HandleFunc(gfa_endpoint, gfa_server_query);
	//http.Handle("/", http.FileServer(http.Dir("js/")));
	http.Handle("/js/", http.StripPrefix("/js/", http.FileServer(http.Dir(gfa_js_dir))));
	if gfa_html_dir != "" {
		http.Handle("/", http.FileServer(http.Dir(gfa_html_dir)));
	}
	fmt.Fprintf(os.Stderr, "[%d] server started at %s\n", time.Now().UnixNano(), gfa_endpoint);

	if gfa_server_socket != "" {
		socket, err := net.Listen("unix", gfa_server_socket);
		if err != nil {
			fmt.Fprintf(os.Stderr, "[%d] could not create socket '%s'\n", time.Now().UnixNano(), gfa_server_socket);
			os.Exit(1);
		}

		c := make(chan os.Signal, 1);
		signal.Notify(c, os.Interrupt, syscall.SIGTERM);
		go func() {
			<-c;
			os.Remove(gfa_server_socket);
			os.Exit(1);
		}()

		http.Serve(socket, nil);
		
	} else {
		http.ListenAndServe(fmt.Sprintf(":%s", gfa_server_port), nil);
	}
}
