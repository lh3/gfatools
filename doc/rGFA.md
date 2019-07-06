## The Graphical Fragment Assembly (GFA) Format

Initially designed and widely used for sequence assembly graphs, GFA can also
encode arbitrary bidirected sequence graphs. GFA primarily consists of two
types of lines: segment lines (S-lines) and link lines (L-lines). An S-line
represents a sequence; an L-line represents a connection between two oriented
segments. The core of GFA can be loosly described by the following grammar:

```txt
<gfa>    <- <line>+
<line>   <- <S-line> | <L-line>
<S-line> <- 'S' <segId> <seq> <tag>*
<L-line> <- 'L' <segId> <strand> <segId> <strand> <cigar> <tag>*
<strand> <- '+' | '-'
```

where fields on each line are TAB delimited; `<cigar>` and `<tag>` follow the
same format as in [SAM][sam]. For details, please see the [GFA spec][gfa1].

## The Reference GFA (rGFA) Format

In GFA, each base can be indexed by a segment ID and an offset on the segment.
This gives us the *segment coordinate* of the base. The segment coordinate is
unstable -- when we split a segment in half, the coordinate is changed. As a
pan-genome graph demands long-term stability, an ordinary GFA is not a good
fit. rGFA address this issue.

rGFA is a strict subset of GFA. It disallows overlaps between segments and
requires three additional tags on each segment. These tags trace the origin of
the segment:

|Tag |Type|Description|
|:--:|:--:|:----------|
|`SN`|`Z` |name of stable sequence from which the segment is derived|
|`SO`|`i` |offset on the stable sequence|
|`SR`|`i` | `0` if the segment is on a linear reference genome; `>0` otherwise|

<img align="right" width="250" src="example1.png"/>

In rGFA, each base in the graph is uniquely indexed by the stable sequence
name and the offset on the stable sequence. This is called the *stable
coordinate* of the base. The stable coordinate never changes as long as bases
remain in the graph.

The figure on the right shows an example rGFA. We can pinpoint a position
such as `chr1:9` in the graph and maps existing annotations onto it. We can
similarly denote a walk or path in the stable coordinate. For example, path
`v1->v2->v3->v4` corresponds to `chr1:0-17` and path `v1->v2->v5->v6`
corresponds to `chr1:0-8=>foo:8-16`. Conversely, an interval on a linear
reference is uniquely represented by a path in the graph. This way rGFA
establishes a stable connection between linear sequences and sequence graphs.

[sam]: https://en.wikipedia.org/wiki/SAM_(file_format)
[gfa1]: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
