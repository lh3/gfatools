function gfa_plot_arrow(ctx, x, y, len, w, color, rev, text)
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
	ctx.lineWidth = 1;
	ctx.stroke();
}

function gfa_plot_conf()
{
	return {
		min_len: 10,
		scale:   10,
		xskip:   15,
		yskip:   30
	};
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
		if (g.idx[v^1][1] == 0 && v0_src < 0)
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
			pred[sub.a[sub.v[i].off + j][0]].push(i);

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
				var end = pos_ij.start + pos_ij.len;
				if (end > max_end) max_end = end;
				pl[pos_ij.level].i = pred[i][j];
				++pl[pos_ij.level].cnt;
			}
			pos[i].start = max_end + conf.xskip;
			// look for an existing level
			var l;
			for (l = 0; l < level_max.length; ++l)
				if (pl[l].cnt < 2 && level_max[l] + conf.xskip <= pos[i].start)
					break;
			pos[i].level = l;
		}
		level_max[pos[i].level] = pos[i].start + pos[i].len;
	}
	var lines = [];
	for (var i = 0; i < pos.length; ++i)
		lines.push([g.seg[sub.v[i].v>>1].name, pos[i].level, pos[i].start, pos[i].len].join("\t"));
//	alert(lines.join("\n"));
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

	for (var i = 0; i < pos.length; ++i) {
		for (var j = 0; j < sub.v[i].n; ++j) {
			ctx.moveTo(pos[i].cx_en, pos[i].cy);
			var k = sub.a[sub.v[i].off + j][0];
			ctx.lineTo(pos[k].cx_st, pos[k].cy);
			ctx.strokeStyle = "#A0A0A0";
			ctx.stroke();
		}
	}
	for (var i = 0; i < pos.length; ++i) {
		var s = g.seg[sub.v[i].v>>1];
		var color = s.rank == 0? '#FF0000' : '#000000';
		var label;
		if (s.rank == 0) label = s.sname + ":" + s.soff + ":" + s.len;
		else label = s.len;
		label = s.len;
		gfa_plot_arrow(ctx, pos[i].cx_st, pos[i].cy, pos[i].len, 4, color, sub.v[i].v&1, label);
	}
}

function gfa_plot(canvas, gfa_text, info_elem)
{
	var g = gfa_parse(gfa_text);
	if (g == null || g.seg.length == 0) {
		info_elem.innerHTML = "Error: failed to parse GFA";
		return;
	}
	var conf = gfa_plot_conf();
	gfa_plot_draw(canvas, conf, g);
}
