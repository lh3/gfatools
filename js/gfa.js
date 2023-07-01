function gfa_seg_add(g, name)
{
	if (name in g.segname) {
		return g.segname[name];
	} else {
		var sid = g.seg.length;
		g.segname[name] = sid;
		g.seg.push({ name:name, len:-1, sname:null, soff:-1, rank:-1 });
		return sid;
	}
}

function gfa_index(g)
{
	var n_vtx = g.seg.length * 2;
	for (var v = 0; v < n_vtx; ++v)
		g.idx[v] = { o:0, n:0 };
	g.arc = g.arc.sort(function(a,b) { return a.v - b.v });
	var st = 0;
	for (var i = 1; i <= g.arc.length; ++i)
		if (i == g.arc.length || g.arc[i].v != g.arc[st].v)
			g.idx[g.arc[st].v] = { o:st, n:i-st }, st = i;
	for (var v = 0; v < n_vtx; ++v) {
		var ov = g.idx[v].o;
		var nv = g.idx[v].n;
		var i0 = -1, n0 = 0;
		for (var i = 0; i < nv; ++i)
			if (g.arc[ov + i].rank == 0)
				++n0, i0 = i;
		if (n0 > 1) g.err |= 2;
		if (i0 > 0) { // then swap [0] and [i0]
			var tmp = g.arc[ov];
			g.arc[ov] = g.arc[ov + i0];
			g.arc[ov + i0] = tmp;
		}
	}
}

function gfa_parse(str)
{
	var g = { seg:[], arc:[], segname:{}, idx:[], walk:[], err:0 };
	var lines = str.split("\n");
	var re_cigar = /(\d+)([MIDSN])/g;
	var re_walk = /([><])([^\s><]+)/g;
	for (var i = 0; i < lines.length; ++i) {
		if (lines[i].length < 5) continue;
		var m, t = lines[i].split("\t");
		if (t[0] == "S") {
			var sid = gfa_seg_add(g, t[1]);
			var s = g.seg[sid];
			if (t[2] != "*") s.len = t[2].length;
			for (var j = 3; j < t.length; ++j) {
				if ((m = /^(LN:i|SN:Z|SO:i|SR:i):(\S+)/.exec(t[j])) == null) continue;
				if (m[1] == "LN:i") s.len = parseInt(m[2]);
				else if (m[1] == "SN:Z") s.sname = m[2];
				else if (m[1] == "SO:i") s.soff = parseInt(m[2]);
				else if (m[1] == "SR:i") s.rank = parseInt(m[2]);
			}
		} else if (t[0] == "L") {
			if (t.length < 5) continue;
			if (t[2] != '+' && t[2] != '-') continue;
			if (t[4] != '+' && t[4] != '-') continue;
			var sid1 = gfa_seg_add(g, t[1]);
			var sid2 = gfa_seg_add(g, t[3]);
			var v = sid1<<1 | (t[2] == '+'? 0 : 1);
			var w = sid2<<1 | (t[4] == '+'? 0 : 1);
			var ov = 0, ow = 0, rank = -1;
			for (var j = 6; j < t.length; ++j)
				if ((m = /^(SR:i):(\S+)/.exec(t[j])) != null)
					rank = parseInt(m[2]);
			if (t.length >= 6) {
				while ((m = re_cigar.exec(t[5])) != null) {
					if (m[2] == 'M' || m[2] == 'D' || m[2] == 'N') ov += parseInt(m[1]);
					if (m[2] == 'M' || m[2] == 'I' || m[2] == 'S') ow += parseInt(m[1]);
				}
			}
			g.arc.push({ v:v, w:w, ov:ov, ow:ow, rank:rank, ori:true });
			g.arc.push({ v:w^1, w:v^1, ov:ow, ow:ov, rank:rank, ori:false });
		} else if (t[0] == "W") {
			if (t.length < 7) continue;
			var walk = { sample:t[1], hap:parseInt(t[2]), sname:t[3], st:-1, en:-1, v:[] };
			if (t[4] != "*") walk.st = parseInt(t[4]);
			if (t[5] != "*") walk.st = parseInt(t[5]);
			while ((m = re_walk.exec(t[6])) != null) {
				if (g.segname[m[2]] != null) {
					var sid = g.segname[m[2]];
					var v = sid<<1 | (m[1] == '>'? 0 : 1);
					walk.v.push(v);
				}
			}
			g.walk.push(walk);
		}
	}

	for (var i = 0; i < g.seg.length; ++i)
		if (g.seg[i].len < 0) g.err |= 1;
	gfa_index(g);
	return g;
}

function gfa_scc1_aux(g)
{
	var n_vtx = g.seg.length * 2;
	var aux = { a:[], index:0 };
	for (var i = 0; i < n_vtx; ++i)
		aux.a.push({ index:-1, start:-1, low:0, i:-1, stack:false });
	return aux;
}

function gfa_scc1(g, aux, v0)
{
	var sub = { v:[], a:[] };
	var ds = [], ts = [];
	ds.push([v0, 0]);
	while (ds.length > 0) {
		var x = ds.pop();
		var v = x[0], i = x[1];
		if (i == 0) { // i is the number of outgoing edges already visited
			aux.a[v].low = aux.a[v].index = aux.index++;
			aux.a[v].stack = true;
			ts.push(v);
		}
		var nv = g.idx[v].n;
		if (i == nv) { // done with v
			if (aux.a[v].low == aux.a[v].index) {
				while (ts.length > 0) {
					var w = ts.pop();
					sub.v.push({ v:w, off:0, n:0 });
					aux.a[w].stack = false;
					if (w == v) break;
				}
			}
			if (ds.length > 0) { // if the DFS stack is not empty, update the top element
				var w = v;
				v = ds[ds.length - 1][0];
				aux.a[v].low = aux.a[v].low < aux.a[w].low? aux.a[v].low : aux.a[w].low;
			}
		} else { // process v's neighbor av[i].w
			var w = g.arc[g.idx[v].o + i].w;
			ds.push([v, i + 1]);
			if (aux.a[w].index == -1 && aux.a[w^1].stack == false)
				ds.push([w, 0]);
			else if (aux.a[w].stack)
				aux.a[v].low = aux.a[v].low < aux.a[w].index? aux.a[v].low : aux.a[w].index;
		}
	}
	// reverse the sub.v[] array
	for (var k = 0; k < sub.v.length>>1; ++k) {
		var x = sub.v[k];
		sub.v[k] = sub.v[sub.v.length - 1 - k];
		sub.v[sub.v.length - 1 - k] = x;
	}
	// fill other fields in sub
	for (var k = 0; k < sub.v.length; ++k)
		aux.a[sub.v[k].v].start = v0, aux.a[sub.v[k].v].i = k;
	for (var k = 0; k < sub.v.length; ++k) {
		var o0 = sub.a.length, v = sub.v[k].v, nv = g.idx[v].n, ov = g.idx[v].o;
		for (var i = 0; i < nv; ++i) {
			var a = g.arc[ov + i];
			if (aux.a[a.w].start == v0)
				sub.a.push({ i:aux.a[a.w].i, arc_off:ov+i, rank:a.rank });
		}
		sub.v[k].off = o0;
		sub.v[k].n = sub.a.length - o0;
		if (sub.v[k].n > 1) { // sort
			var b = sub.a.slice(o0).sort(function(x,y) { return x[0] - y[0] });
			for (var i = 0; i < b.length; ++i)
				sub.a[o0 + i] = b[i];
		}
	}
	return sub;
}

function gfa_scc1_string(g, sub)
{
	var lines = [];
	for (var i = 0; i < sub.v.length; ++i) {
		var v = sub.v[i].v;
		var t = ["[" + i + "]", v];
		t.push("><"[v&1] + g.seg[v>>1].name, sub.v[i].n);
		if (sub.v[i].n > 0) {
			var s = [];
			for (var j = 0; j < sub.v[i].n; ++j)
				s.push(sub.a[sub.v[i].off + j].i);
			t.push(s.join(","));
		}
		lines.push(t.join("\t"));
	}
	return lines.join("\n");
}

function gfa_string(g)
{
	var lines = [];
	for (var i = 0; i < g.seg.length; ++i) {
		var s = g.seg[i], t = ['S', s.name, '*'];
		if (s.len >= 0) t.push('LN:i:' + s.len);
		if (s.sname != null && s.soff >= 0)
			t.push('SN:Z:' + s.sname, 'SO:i:' + s.soff);
		if (s.rank >= 0) t.push('SR:i:' + s.rank);
		lines.push(t.join("\t"));
	}
	for (var i = 0; i < g.arc.length; ++i) {
		var a = g.arc[i];
		if (!a.ori) continue;
		var t = ['L', g.seg[a.v>>1].name, a.v&1? '-' : '+', g.seg[a.w>>1].name, a.w&1? '-' : '+'];
		if (a.ov == 0 && a.ow == 0) t.push('0M');
		else t.push(a.ov + ':' + a.ow);
		if (a.rank >= 0) t.push('SR:i:' + a.rank);
		lines.push(t.join("\t"));
	}
	return lines.join("\n");
}
