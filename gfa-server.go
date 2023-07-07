package main

import (
	"os"
	"fmt"
	"unsafe"
	"strconv"
	"net/http"
	"strings"
	"time"
)

/*
#cgo LDFLAGS: kalloc.o gfa-base.o gfa-io.o gfa-util.o -lz -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfa-priv.h"

char *gfa_extract(gfa_t *g, const char *str, int step)
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
	gfa_walk_flip(f);
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
var gfa_graph *C.gfa_t;

func gfa_server_query(w http.ResponseWriter, r *http.Request) {
	r.ParseForm();
	start_time := time.Now().UnixNano();
	fmt.Fprintf(os.Stderr, "[%d] got request: %s\n", start_time, r.Form);
	defer fmt.Fprintf(os.Stderr, "[%d] responded %d\n", time.Now().UnixNano(), start_time);
	if len(r.Form) == 0 {
		fmt.Fprintln(w, "Hello!");
		return;
	}
	step := 0;
	if len(r.Form["step"]) > 0 { // set radius
		i, _ := strconv.Atoi(r.Form["step"][0]);
		if i < 0 {
			http.Error(w, "400 Bad Request: 'r' shouldn't be negative", 400);
			return;
		}
		step = i;
	}
	if len(r.Form["q"]) > 0 {
		cstr := C.CString(r.Form["q"][0]);
		out := C.gfa_extract(gfa_graph, cstr, C.int(step));
		C.free(unsafe.Pointer(cstr));
		ret := C.GoString(out);
		C.free(unsafe.Pointer(out));
		if len(r.Form["plot"]) > 0 {
			w.Header().Set("Content-Type", "text/html; charset=utf-8");
			fmt.Fprintln(w, "<title>GFA view</title>");
			fmt.Fprintln(w, "<style type='text/css'>#canvas_graph,#canvas_walk { border: 1px solid #000; }</style>");
			fmt.Fprintln(w, "<script language='JavaScript' src='gfa.js'></script>");
			fmt.Fprintln(w, "<script language='JavaScript' src='gfa-plot.js'></script>");
			fmt.Fprintln(w, "<body>");
			fmt.Fprintln(w, "<p><canvas id='canvas_graph' width='800' height='100'></canvas></p>");
			fmt.Fprintln(w, "<p><canvas id='canvas_walk' width='800' height='100'></canvas></p>");
			fmt.Fprintln(w, "<textarea id='gfa-text' readonly rows='15' cols='110'>");
			fmt.Fprintf(w, ret);
			fmt.Fprintln(w, "</textarea>");
			fmt.Fprintln(w, "<script language='JavaScript'>");
			fmt.Fprintln(w, "var gfa_conf = gfa_plot_conf();");
			fmt.Fprintln(w, "var gfa_text = document.getElementById('gfa-text').value;");
			fmt.Fprintln(w, "var gfa_graph = gfa_parse(gfa_text);");
			fmt.Fprintln(w, "gfa_conf.merged = true; gfa_conf.uniq_walk = true;");
			fmt.Fprintln(w, "gfa_plot_graph(document.getElementById('canvas_graph'), gfa_conf, gfa_graph);");
			fmt.Fprintln(w, "gfa_plot_walk(document.getElementById('canvas_walk'), gfa_conf, gfa_graph);");
			fmt.Fprintln(w, "</script>");
			fmt.Fprintln(w, "</body>");
		} else {
			fmt.Fprintf(w, "%s", ret);
		}
	}
}

/*****************
 * Main function *
 *****************/

func main() {
	if os.Getenv("PORT") != "" {
		gfa_server_port = os.Getenv("PORT");
	}
	// parse command line options
	for {
		opt, arg := getopt(os.Args, "p:");
		if opt == 'p' {
			gfa_server_port = arg;
		} else if opt < 0 {
			break;
		}
	}
	if optind == len(os.Args) {
		fmt.Fprintln(os.Stderr, "Usage: gfa-server [options] <graph.gfa>");
		fmt.Fprintln(os.Stderr, "Options:");
		fmt.Fprintf(os.Stderr, "  -p INT    port number [%s or from $PORT env]\n", gfa_server_port);
		os.Exit(1);
	}

	fn := C.CString(os.Args[optind]);
	defer C.free(unsafe.Pointer(fn));
	gfa_graph = C.gfa_read(fn);
	defer C.gfa_destroy(gfa_graph);

	root := "/gfasub";
	http.HandleFunc(root, gfa_server_query);
	http.Handle("/", http.FileServer(http.Dir("js/")));
	//http.Handle("/static/", http.StripPrefix("/static/", http.FileServer(http.Dir("js/"))));
	fmt.Fprintf(os.Stderr, "[%d] server started at %s\n", time.Now().UnixNano(), root);
	http.ListenAndServe(fmt.Sprintf(":%s", gfa_server_port), nil);
}
