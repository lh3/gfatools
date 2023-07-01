function gfa_plot_arrow(ctx, x, y, len, w, color, rev, text, lw)
{
	ctx.font = "10px Arial";
	if (text != null) {
		ctx.strokeStyle = "#000000";
		ctx.textAlign = "center";
		ctx.fillText(text, x + len/2, y - w - 2);
	}
	ctx.beginPath();
	ctx.moveTo(x, y);
	if (rev == null || !rev) {
		ctx.lineTo(x - w, y - w);
		ctx.lineTo(x - w + len, y - w);
		ctx.lineTo(x + len, y);
		ctx.lineTo(x - w + len, y + w);
		ctx.lineTo(x - w, y + w);
	} else {
		ctx.lineTo(x + w, y - w);
		ctx.lineTo(x + w + len, y - w);
		ctx.lineTo(x + len, y);
		ctx.lineTo(x + w + len, y + w);
		ctx.lineTo(x + w, y + w);
	}
	ctx.lineTo(x, y);
	ctx.strokeStyle = color;
	ctx.lineWidth = lw == null? 1 : lw;
	ctx.stroke();
}

function gfa_plot_conf()
{
	return {
		label:   "name", // or "len"
		min_len: 10,
		scale:   10,
		xskip:   15,
		yskip:   30
	};
}

var gfa_conf = gfa_plot_conf();

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

function gfa_plot_cal_pos(conf, g)
{
	function cal_length(len, min_len, scale) {
		return Math.floor(min_len + Math.log(len + 1) / Math.log(10) * scale + .499);
	}

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
		pos[i].len = cal_length(g.seg[sub.v[i].v>>1].len, conf.min_len, conf.scale);
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

function gfa_plot_draw(canvas, conf, g)
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

	// draw nodes
	ctx.globalAlpha = 1.0;
	for (var i = 0; i < pos.length; ++i) {
		var s = g.seg[sub.v[i].v>>1];
		var label;
		if (conf.label == "name") label = s.name;
		else if (conf.label == "length") label = s.len;
		var color = r2c[s.rank] != null? r2c[s.rank] : "#000000";
		gfa_plot_arrow(ctx, pos[i].cx_st, pos[i].cy, pos[i].len, 4, color, sub.v[i].v&1, label, s.rank == 0? 1.5 : 1);
	}
}

function gfa_plot(canvas, gfa_text, info_elem)
{
	var g = gfa_parse(gfa_text);
	if (g == null || g.seg.length == 0) {
		info_elem.innerHTML = "Error: failed to parse GFA";
		return;
	}
	gfa_plot_draw(canvas, gfa_conf, g);
}
