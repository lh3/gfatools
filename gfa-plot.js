function gfa_plot_conf()
{
	return {
		label:      "name", // or "len"
		merge_walk: true,
		uniq_walk:  true,
		font_size:  9,
		min_len:    10,
		scale:      10,
		h_arrow:    4,
		xskip:      15,
		yskip:      30
	};
}

var gfa_conf = gfa_plot_conf();
var gfa_obj = null;

function gfa_plot_arrow(ctx, x, y, len, w, rev, text, fs, color_stroke, color_fill, lw)
{
	ctx.font = fs? fs + "px mono" : "9px mono";
	var reg = new Path2D();
	reg.moveTo(x, y);
	if (rev == null || !rev) {
		reg.lineTo(x - w, y - w);
		reg.lineTo(x - w + len, y - w);
		reg.lineTo(x + len, y);
		reg.lineTo(x - w + len, y + w);
		reg.lineTo(x - w, y + w);
	} else {
		reg.lineTo(x + w, y - w);
		reg.lineTo(x + w + len, y - w);
		reg.lineTo(x + len, y);
		reg.lineTo(x + w + len, y + w);
		reg.lineTo(x + w, y + w);
	}
	reg.lineTo(x, y);
	reg.closePath();
	if (color_fill) {
		ctx.fillStyle = color_fill;
		ctx.fill(reg);
	}
	if (color_stroke) {
		ctx.strokeStyle = color_stroke;
		ctx.stroke(reg);
	}
	if (text != null) {
		ctx.textAlign = "center";
		ctx.fillText(text, x + len/2, y - w - 2);
		ctx.fillStyle = color_fill? color_fill : "#000000";
		ctx.stroke();
	}
}

function gfa_plot_rank2color(g)
{
	var color = [ "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf" ];
	var rank = [];
	for (var i = 0; i < g.seg.length; ++i)
		if (g.seg[i].rank >= 0) rank.push(g.seg[i].rank);
	for (var i = 0; i < g.arc.length; ++i)
		if (g.arc[i].rank >= 0) rank.push(g.arc[i].rank);
	if (rank.length == 0) return [];
	rank = rank.sort();
	var shrink = [rank[0]];
	for (var i = 1; i < rank.length; ++i)
		if (rank[i] != rank[i-1])
			shrink.push(rank[i]);
	var r2c = {};
	for (var i = 0; i < shrink.length && i < color.length; ++i)
		r2c[shrink[i]] = color[i];
	return r2c;
}

function gfa_plot_find_v0(g) // FIXME: not general
{
	var n_vtx = g.seg.length * 2, v0_ref = -1, v0_src = -1;
	for (var v = 0; v < n_vtx; ++v) {
		if (g.seg[v>>1].rank == 0) {
			if (g.seg[v>>1].snid < 0 || g.seg[v>>1].soff < 0) {
				info_elem.innerHTML = "Error: a rank-0 segment without stable coordinates";
				return;
			}
			if (v0_ref < 0 || g.seg[v>>1].soff < g.seg[v0_ref>>1].soff)
				v0_ref = v;
		}
		if (g.idx[v^1].n == 0 && v0_src < 0)
			v0_src = v;
	}
	return v0_ref >= 0? v0_ref : v0_src >= 0? v0_src : 0;
}

function gfa_plot_cal_length(len, min_len, scale)
{
	return Math.floor(min_len + Math.log(len + 1) / Math.log(10) * scale + .499);
}

function gfa_plot_cal_pos(conf, g)
{
	var v0 = gfa_plot_find_v0(g);
	var aux = gfa_scc1_aux(g);
	var sub = gfa_scc1(g, aux, v0);
	var pred = [];
	for (var i = 0; i < sub.v.length; ++i) pred[i] = [];
	for (var i = 0; i < sub.v.length; ++i)
		for (var j = 0; j < sub.v[i].n; ++j)
			pred[sub.a[sub.v[i].off + j].i].push(i);

	var l = 0, level_max = [], pos = [];
	for (var i = 0; i < sub.v.length; ++i) {
		pos[i] = { level:-1, start:-1, len:0 };
		pos[i].len = gfa_plot_cal_length(g.seg[sub.v[i].v>>1].len, conf.min_len, conf.scale);
	}
	for (var i = 0; i < sub.v.length; ++i) {
		var t = sub.v[i];
		if (pred[i].length == 0) { // a source node
			pos[i].start = 0;
			pos[i].level = level_max.length;
		} else {
			var max_end = -1, pl = [];
			for (var j = 0; j < level_max.length; ++j) pl[j] = { cnt:0, i:-1 };
			for (var j = 0; j < pred[i].length; ++j) {
				var pos_ij = pos[pred[i][j]];
				if (pos_ij.level >= 0) {
					var end = pos_ij.start + pos_ij.len;
					if (end > max_end) max_end = end;
					pl[pos_ij.level].i = pred[i][j];
					pl[pos_ij.level].end = pos[pred[i][j]].start + pos[pred[i][j]].len;
					++pl[pos_ij.level].cnt;
				}
			}
			pos[i].start = max_end + conf.xskip;
			// look for an existing level
			var l;
			for (l = 0; l < level_max.length; ++l) {
				if (pl[l].cnt > 1 || level_max[l] + conf.xskip > pos[i].start) continue;
				if (pl[l].cnt == 0) break;
				if (pl[l].end == level_max[l]) break;
			}
			pos[i].level = l;
		}
		level_max[pos[i].level] = pos[i].start + pos[i].len;
	}
	return [sub, pos];
}

function gfa_plot_graph(canvas, conf, g)
{
	var ctx = canvas.getContext("2d");
	ctx.clearRect(0, 0, canvas.width, canvas.height);

	var ret = gfa_plot_cal_pos(conf, g);
	var sub = ret[0], pos = ret[1];

	for (var i = 0; i < pos.length; ++i) {
		pos[i].cx_st = conf.xskip + pos[i].start;
		pos[i].cx_en = pos[i].cx_st + pos[i].len;
		pos[i].cy = conf.yskip * (pos[i].level + 1);
	}

	var max_w = 0, max_h = 0;
	for (var i = 0; i < pos.length; ++i) {
		var x = pos[i].cx_en + conf.xskip;
		max_w = max_w > x? max_w : x;
		var y = pos[i].cy + conf.yskip;
		max_h = max_h > y? max_h : y;
	}
	canvas.width = max_w, canvas.height = max_h;

	ctx = canvas.getContext("2d");
	ctx.translate(0.5, 0.5);

	var r2c = gfa_plot_rank2color(g);

	// draw edges
	ctx.globalAlpha = 0.8;
	for (var i = 0; i < pos.length; ++i) {
		for (var j = 0; j < sub.v[i].n; ++j) {
			ctx.beginPath();
			ctx.moveTo(pos[i].cx_en, pos[i].cy);
			var k = sub.a[sub.v[i].off + j].i;
			ctx.lineTo(pos[k].cx_st, pos[k].cy);
			var r = sub.a[sub.v[i].off + j].rank;
			ctx.strokeStyle = r >= 0 && r2c[r] != null? r2c[r] : "#A0A0A0";
			ctx.stroke();
		}
	}
	ctx.globalAlpha = 1.0;

	// draw nodes
	ctx.globalAlpha = 1.0;
	for (var i = 0; i < pos.length; ++i) {
		var s = g.seg[sub.v[i].v>>1];
		var label;
		if (conf.label == "name") label = s.name;
		else if (conf.label == "length") label = s.len;
		var color_stroke = s.rank >= 0 && r2c[s.rank] != null? r2c[s.rank] : null;
		var lw = s.rank == 0? 1.5 : 1;
		gfa_plot_arrow(ctx, pos[i].cx_st, pos[i].cy, pos[i].len, conf.h_arrow, sub.v[i].v&1, label, conf.font_size, color_stroke, s.color, lw);
	}
}

/*
 * Plot walks
 */
function gfa_int_hash(x)
{
	x = ((x >> 16) ^ x) * 0x45d9f3b & 0xffffffff;
	x = ((x >> 16) ^ x) * 0x45d9f3b & 0xffffffff;
	return (x >> 16) ^ x;
}

function gfa_walk_gen(g, merge, uniq) // filter or combine walks
{
	var walk = [], tmp = [];
	if (uniq) { // only choose samples with one walk
		var t2 = g.walk.sort(function(x,y) { return x.asm == y.asm? 0 : x.asm < y.asm? -1 : 1 });
		var i, i0;
		for (i0 = 0, i = 1; i <= t2.length; ++i) {
			if (i == t2.length || t2[i].asm != t2[i0].asm) {
				if (i - i0 == 1) tmp.push(t2[i0]);
				i0 = i;
			}
		}
	} else { // choose all walks
		for (var i = 0; i < g.walk.length; ++i)
			tmp.push(g.walk[i]);
	}
	for (var i = 0; i < tmp.length; ++i) { // compute hash
		var w = tmp[i];
		var ww = { label:w.asm, hash:0, n:1, v:w.v };
		var hash = 0;
		for (var j = 0; j < w.v.length; ++j)
			hash = (hash + gfa_int_hash(w.v[j])) & 0xffffffff;
		ww.hash = w.v.length << 32 | hash;
		walk.push(ww);
	}
	if (merge) { // merge identical walks; FIXME: this only checks hash
		walk = walk.sort(function(x,y) { return x.hash - y.hash; });
		var i0, i, k;
		for (i0 = 0, i = 1, k = 0; i <= walk.length; ++i) {
			if (i == walk.length || walk[i0].hash != walk[i].hash) {
				walk[k] = walk[i0];
				walk[k++].n = i - i0;
				i0 = i;
			}
		}
		walk.length = k;
		walk = walk.sort(function(x,y) { return y.n - x.n; });
		for (i = 0; i < walk.length; ++i) // reassign sample name
			walk[i].label = "" + walk[i].n;
	}
	return walk;
}

function gfa_plot_walk(canvas, conf, g)
{
	if (g.walk.length == 0) return;

	var seg_aux = [];
	for (var i = 0; i < g.seg.length; ++i) {
		seg_aux[i] = {};
		seg_aux[i].clen = gfa_plot_cal_length(g.seg[i].len, conf.min_len, conf.scale);
	}

	var walk = gfa_walk_gen(g, conf.merge_walk, conf.uniq_walk);

	var max_len = 0, max_label_len = 0;
	for (var i = 0; i < walk.length; ++i) {
		var w = walk[i], len = 0;
		max_label_len = max_label_len > w.label.length? max_label_len : w.label.length;
		for (var j = 0; j < w.v.length; ++j)
			len += seg_aux[w.v[j]>>1].clen;
		len += (w.v.length - 1) * conf.xskip;
		max_len = max_len > len? max_len : len;
	}
	var off_x = conf.xskip + (max_label_len + 1) * conf.font_size + conf.xskip;
	var max_w = off_x + max_len + conf.xskip;
	var max_h = (walk.length + 1) * conf.yskip;
	canvas.width = max_w, canvas.height = max_h;

	var ctx = canvas.getContext("2d");

	// draw segments
	var cy = conf.yskip;
	for (var i = 0; i < walk.length; ++i) {
		var w = walk[i], cx = off_x, cx_pre;
		for (var j = 0; j < w.v.length; ++j) {
			var label, sid = w.v[j]>>1;
			var s = g.seg[sid];
			if (conf.label == "name") label = s.name;
			else if (conf.label == "length") label = s.len;
			gfa_plot_arrow(ctx, cx, cy, seg_aux[sid].clen, conf.h_arrow, w.v[j]&1, label, conf.font_size, null, s.color, 1);
			cx_pre = cx + seg_aux[sid].clen;
			cx += seg_aux[sid].clen + conf.xskip;
		}
		cy += conf.yskip;
	}

	// draw lines
	cy = conf.yskip;
	ctx.lineWidth = 0.2;
	ctx.strokeStyle = "#404040";
	for (var i = 0; i < walk.length; ++i) {
		var w = walk[i], cx = off_x, cx_pre = 0;
		for (var j = 0; j < w.v.length; ++j) {
			var sid = w.v[j]>>1;
			if (j > 0) {
				ctx.beginPath();
				ctx.moveTo(cx_pre, cy);
				ctx.lineTo(cx, cy);
				ctx.stroke();
			}
			cx_pre = cx + seg_aux[sid].clen;
			cx += seg_aux[sid].clen + conf.xskip;
		}
		cy += conf.yskip;
	}

	// draw labels
	cy = conf.yskip + (conf.font_size>>1);
	for (var i = 0; i < walk.length; ++i) {
		var w = walk[i], cx = off_x, cx_pre = 0;
		ctx.fillStyle = "#000000";
		ctx.textAlign = "left";
		ctx.fillText(w.label, conf.xskip, cy);
		cy += conf.yskip;
	}
}
