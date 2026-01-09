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
		yskip:      30,
		// options for grouping
		group_data:    null,   // parsed TSV: { sample -> {group, group_order, sample_order} }
		group_mode:    false,  // enable group panel mode
		group_gap:     20      // vertical gap between group panels
	};
}

function gfa_parse_group_tsv(text)
{
	var result = {};
	var lines = text.trim().split("\n");
	if (lines.length < 2) return result; // need header + at least 1 data row

	// parse header to find column indices
	var header = lines[0].split("\t");
	var col_sample = -1, col_group = -1, col_group_order = -1, col_sample_order = -1;
	for (var i = 0; i < header.length; ++i) {
		var h = header[i].toLowerCase().trim();
		if (h == "sample") col_sample = i;
		else if (h == "group") col_group = i;
		else if (h == "group_order") col_group_order = i;
		else if (h == "sample_order") col_sample_order = i;
	}

	if (col_sample < 0 || col_group < 0) {
		console.error("TSV must have 'sample' and 'group' columns");
		return null;
	}

	// parse data rows
	for (var i = 1; i < lines.length; ++i) {
		if (lines[i].trim() == "") continue;
		var t = lines[i].split("\t");
		var sample = t[col_sample] ? t[col_sample].trim() : "";
		var group = t[col_group] ? t[col_group].trim() : "";
		if (sample == "" || group == "") continue;

		var entry = { group: group };
		if (col_group_order >= 0 && t[col_group_order])
			entry.group_order = parseInt(t[col_group_order].trim());
		if (col_sample_order >= 0 && t[col_sample_order])
			entry.sample_order = parseInt(t[col_sample_order].trim());

		result[sample] = entry;
	}

	return result;
}

function gfa_validate_groups(g, group_data)
{
	if (group_data == null) return { valid: false, error: "No group data provided" };

	var missing = [];
	var seen_samples = {};

	// collect unique sample names from walks
	for (var i = 0; i < g.walk.length; ++i) {
		var sample = g.walk[i].sample;
		if (seen_samples[sample]) continue;
		seen_samples[sample] = true;

		if (!(sample in group_data)) {
			missing.push(sample);
		}
	}

	if (missing.length > 0) {
		return {
			valid: false,
			error: "Samples missing from TSV: " + missing.join(", "),
			missing: missing
		};
	}

	return { valid: true };
}

function gfa_canvas_hi_res(canvas, width, height)
{
    ratio = window.devicePixelRatio;
    canvas.width = width * ratio;
    canvas.height = height * ratio;
    canvas.style.width = width + "px";
    canvas.style.height = height + "px";
    canvas.getContext("2d").scale(ratio, ratio);
}

var gfa_conf = gfa_plot_conf();
var gfa_obj = null;

function gfa_plot_arrow(ctx, x, y, len, w, rev, text, fs, color_stroke, color_fill)
{
	ctx.font = fs? fs + "px monospace" : "9px monospace";
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
		ctx.strokeStyle = ctx.fillStyle = color_stroke;
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
	gfa_canvas_hi_res(canvas, max_w, max_h);

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
		gfa_plot_arrow(ctx, pos[i].cx_st, pos[i].cy, pos[i].len, conf.h_arrow, sub.v[i].v&1, label, conf.font_size, color_stroke, s.color);
	}
}

/*
 * Plot walks
 */
function gfa_int_hash(x)
{
	x = Math.imul((x >> 16) ^ x, 0x45d9f3b) & 0xffffffff;
	x = Math.imul((x >> 16) ^ x, 0x45d9f3b) & 0xffffffff;
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
		var ww = { label:w.asm, hash:0, n:1, v:w.v, lof:w.lof };
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

function gfa_walk_gen_grouped(g, group_data, merge, uniq)
{
	// Step 1: Filter walks by haplotype (only haplotypes with single walk if uniq=true)
	var tmp = [];
	if (uniq) {
		var t2 = g.walk.slice().sort(function(x,y) {
			return x.asm == y.asm ? 0 : x.asm < y.asm ? -1 : 1;
		});
		var i, i0;
		for (i0 = 0, i = 1; i <= t2.length; ++i) {
			if (i == t2.length || t2[i].asm != t2[i0].asm) {
				if (i - i0 == 1) tmp.push(t2[i0]);
				i0 = i;
			}
		}
	} else {
		for (var i = 0; i < g.walk.length; ++i)
			tmp.push(g.walk[i]);
	}

	// Step 2: Group walks by group name
	var groups = {};  // group_name -> { walks: [], group_order: number }
	for (var i = 0; i < tmp.length; ++i) {
		var w = tmp[i];
		var sample = w.sample;
		if (!(sample in group_data)) continue;

		var ginfo = group_data[sample];
		var gname = ginfo.group;

		if (!(gname in groups)) {
			groups[gname] = {
				name: gname,
				group_order: ginfo.group_order != null ? ginfo.group_order : 999999,
				walks: []
			};
		}

		// compute hash for this walk
		var hash = 0;
		for (var j = 0; j < w.v.length; ++j)
			hash = (hash + gfa_int_hash(w.v[j])) & 0xffffffff;
		hash = w.v.length << 32 | hash;

		groups[gname].walks.push({
			sample: sample,
			asm: w.asm,  // sample#haplotype label
			sample_order: ginfo.sample_order != null ? ginfo.sample_order : 999999,
			hash: hash,
			v: w.v,
			lof: w.lof
		});
	}

	// Step 3: Sort groups by group_order, then alphabetically
	var group_list = [];
	for (var gname in groups) group_list.push(groups[gname]);
	group_list.sort(function(a, b) {
		if (a.group_order != b.group_order) return a.group_order - b.group_order;
		return a.name < b.name ? -1 : a.name > b.name ? 1 : 0;
	});

	// Step 4: Within each group, process haplotypes
	for (var gi = 0; gi < group_list.length; ++gi) {
		var grp = group_list[gi];
		grp.sample_count = grp.walks.length;

		if (merge) {
			// Collapse identical haplotypes by hash
			grp.walks.sort(function(a, b) { return a.hash - b.hash; });
			var collapsed = [];
			var i0, i;
			for (i0 = 0, i = 1; i <= grp.walks.length; ++i) {
				if (i == grp.walks.length || grp.walks[i0].hash != grp.walks[i].hash) {
					collapsed.push({
						label: "" + (i - i0),
						count: i - i0,
						v: grp.walks[i0].v,
						lof: grp.walks[i0].lof
					});
					i0 = i;
				}
			}
			// sort by count descending
			collapsed.sort(function(a, b) { return b.count - a.count; });
			grp.haplotypes = collapsed;
		} else {
			// Keep individual haplotypes with sample#haplotype labels
			// sort by sample_order, then alphabetically
			grp.walks.sort(function(a, b) {
				if (a.sample_order != b.sample_order) return a.sample_order - b.sample_order;
				return a.asm < b.asm ? -1 : a.asm > b.asm ? 1 : 0;
			});
			var haplotypes = [];
			for (var i = 0; i < grp.walks.length; ++i) {
				haplotypes.push({
					label: grp.walks[i].asm,
					count: 1,
					v: grp.walks[i].v,
					lof: grp.walks[i].lof
				});
			}
			grp.haplotypes = haplotypes;
		}
		delete grp.walks;
	}

	return group_list;
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
	gfa_canvas_hi_res(canvas, max_w, max_h);

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
			if (w.lof == null || w.lof.length == 0 || w.lof[j] == 0)
				gfa_plot_arrow(ctx, cx, cy, seg_aux[sid].clen, conf.h_arrow, w.v[j]&1, label, conf.font_size, null, s.color);
			else
				gfa_plot_arrow(ctx, cx, cy, seg_aux[sid].clen, conf.h_arrow, w.v[j]&1, label, conf.font_size, s.color, null);
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

function gfa_plot_walk_grouped(canvas, conf, g)
{
	if (g.walk.length == 0) return;
	if (conf.group_data == null) return;

	// validate first
	var validation = gfa_validate_groups(g, conf.group_data);
	if (!validation.valid) {
		console.error(validation.error);
		return;
	}

	// generate grouped data
	var groups = gfa_walk_gen_grouped(g, conf.group_data, conf.merge_walk, conf.uniq_walk);
	if (groups.length == 0) return;

	// precompute segment lengths
	var seg_aux = [];
	for (var i = 0; i < g.seg.length; ++i) {
		seg_aux[i] = {};
		seg_aux[i].clen = gfa_plot_cal_length(g.seg[i].len, conf.min_len, conf.scale);
	}

	// calculate max walk length and max label length across all groups
	var max_len = 0, max_label_len = 0;
	for (var gi = 0; gi < groups.length; ++gi) {
		var grp = groups[gi];
		// header label length
		var header = grp.name + " (n=" + grp.sample_count + ", " + grp.haplotypes.length + " haplotypes)";
		max_label_len = Math.max(max_label_len, header.length);

		for (var i = 0; i < grp.haplotypes.length; ++i) {
			var h = grp.haplotypes[i];
			max_label_len = Math.max(max_label_len, h.label.length);

			var len = 0;
			for (var j = 0; j < h.v.length; ++j)
				len += seg_aux[h.v[j] >> 1].clen;
			len += (h.v.length - 1) * conf.xskip;
			max_len = Math.max(max_len, len);
		}
	}

	// calculate total height
	var total_rows = 0;
	for (var gi = 0; gi < groups.length; ++gi) {
		total_rows += 1;  // header row
		total_rows += groups[gi].haplotypes.length;
	}
	total_rows += groups.length - 1;  // gaps between groups (in yskip units)

	var off_x = conf.xskip + (max_label_len + 1) * (conf.font_size * 0.6) + conf.xskip;
	var max_w = off_x + max_len + conf.xskip;
	var max_h = (total_rows + 1) * conf.yskip + (groups.length - 1) * conf.group_gap;
	gfa_canvas_hi_res(canvas, max_w, max_h);

	var ctx = canvas.getContext("2d");
	var cy = conf.yskip;

	for (var gi = 0; gi < groups.length; ++gi) {
		var grp = groups[gi];

		// Draw group header
		var header = grp.name + " (n=" + grp.sample_count + ", " + grp.haplotypes.length + " haplotypes)";
		ctx.fillStyle = "#000000";
		ctx.font = "bold " + conf.font_size + "px monospace";
		ctx.textAlign = "left";
		ctx.fillText(header, conf.xskip, cy + (conf.font_size >> 1));
		cy += conf.yskip;

		// Draw haplotypes
		ctx.font = conf.font_size + "px monospace";
		for (var i = 0; i < grp.haplotypes.length; ++i) {
			var h = grp.haplotypes[i];
			var cx = off_x, cx_pre;

			// draw segments
			for (var j = 0; j < h.v.length; ++j) {
				var label, sid = h.v[j] >> 1;
				var s = g.seg[sid];
				if (conf.label == "name") label = s.name;
				else if (conf.label == "length") label = s.len;
				if (h.lof == null || h.lof.length == 0 || h.lof[j] == 0)
					gfa_plot_arrow(ctx, cx, cy, seg_aux[sid].clen, conf.h_arrow, h.v[j] & 1, label, conf.font_size, null, s.color);
				else
					gfa_plot_arrow(ctx, cx, cy, seg_aux[sid].clen, conf.h_arrow, h.v[j] & 1, label, conf.font_size, s.color, null);
				cx_pre = cx + seg_aux[sid].clen;
				cx += seg_aux[sid].clen + conf.xskip;
			}

			// draw connecting lines
			cx = off_x;
			ctx.lineWidth = 0.2;
			ctx.strokeStyle = "#404040";
			for (var j = 0; j < h.v.length; ++j) {
				var sid = h.v[j] >> 1;
				if (j > 0) {
					ctx.beginPath();
					ctx.moveTo(cx_pre, cy);
					ctx.lineTo(cx, cy);
					ctx.stroke();
				}
				cx_pre = cx + seg_aux[sid].clen;
				cx += seg_aux[sid].clen + conf.xskip;
			}

			// draw label (count when merged, sample#haplotype when not merged)
			ctx.fillStyle = "#000000";
			ctx.textAlign = "left";
			ctx.fillText(h.label, conf.xskip, cy + (conf.font_size >> 1));

			cy += conf.yskip;
		}

		// add gap before next group
		if (gi < groups.length - 1) {
			cy += conf.group_gap;
		}
	}
}

function gfa_plot_walk_auto(canvas, conf, g)
{
	if (conf.group_mode && conf.group_data != null) {
		gfa_plot_walk_grouped(canvas, conf, g);
	} else {
		gfa_plot_walk(canvas, conf, g);
	}
}

/*
 * Auto-detect and load groups.tsv for grouped visualization
 */
var gfa_groups_cache = null;      // cached group data
var gfa_groups_fetched = false;   // whether we've attempted to fetch

function gfa_fetch_groups(callback)
{
	if (gfa_groups_fetched) {
		callback(gfa_groups_cache);
		return;
	}

	var xhr = new XMLHttpRequest();
	xhr.open("GET", "/groups.tsv", true);
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4) {
			gfa_groups_fetched = true;
			if (xhr.status == 200 && xhr.responseText.length > 0) {
				gfa_groups_cache = gfa_parse_group_tsv(xhr.responseText);
				if (gfa_groups_cache && Object.keys(gfa_groups_cache).length > 0) {
					console.log("groups.tsv loaded: " + Object.keys(gfa_groups_cache).length + " samples");
				} else {
					gfa_groups_cache = null;
				}
			}
			callback(gfa_groups_cache);
		}
	};
	xhr.send();
}

// Inject "Group by population" checkbox into the page
var gfa_group_checkbox_injected = false;
function gfa_inject_group_checkbox() {
	if (gfa_group_checkbox_injected) return;
	gfa_group_checkbox_injected = true;

	// Find the replot button to insert before it
	var replotBtn = document.querySelector('input[value="Replot"]');
	if (replotBtn) {
		var span = document.createElement('span');
		span.innerHTML = '<input type="checkbox" id="group_view"/> group by population \uFF5C ';
		replotBtn.parentNode.insertBefore(span, replotBtn);
	}
}

// Wrap gfa_plot_walk to support group view toggle
var _gfa_plot_walk_original = gfa_plot_walk;
gfa_plot_walk = function(canvas, conf, g) {
	// Inject checkbox on first call
	gfa_inject_group_checkbox();

	// Check if group view is enabled
	var groupCheckbox = document.getElementById('group_view');
	var groupEnabled = groupCheckbox && groupCheckbox.checked;

	if (!groupEnabled) {
		// Use original logic when group view is disabled
		_gfa_plot_walk_original(canvas, conf, g);
		return;
	}

	// Group view enabled - fetch and apply groups
	gfa_fetch_groups(function(group_data) {
		if (group_data != null) {
			var validation = gfa_validate_groups(g, group_data);
			if (validation.valid) {
				conf.group_data = group_data;
				conf.group_mode = true;
				gfa_plot_walk_grouped(canvas, conf, g);
				return;
			} else {
				console.warn("groups.tsv validation failed: " + validation.error);
				// Fall back to original
				_gfa_plot_walk_original(canvas, conf, g);
			}
		} else {
			// No groups file, use original
			_gfa_plot_walk_original(canvas, conf, g);
		}
	});
};
