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

<img align="right" width="220" src="example1.png"/>

In GFA, each base can be indexed by a segment ID and an offset on the segment,
which is called the *segment coordinate*. The segment coordinate is unstable --
when we split a segment in half, the segment coordinate is changed. An ordinary
GFA is not a good format to store a pan-genome graph which requires a stable
coordinate system. The rGFA format addresses this issue.

rGFA is a strict subset of GFA. It requires three additional tags on each
segment:

|Tag |Type|Description                                                       |
|:--:|:--:|:-----------------------------------------------------------------|
|`SN`|`Z` |Name of stable sequence from which the segment is derived         |
|`SO`|`i` |Offset on the stable sequence                                     |
|`SR`|`i` |`0` if the segment is on a linear reference genome; `>0` otherwise|

The figure on the right shows an example rGFA.

[miniasm]: https://github.com/lh3/miniasm/
[canu]: https://github.com/marbl/canu
[spades]: https://github.com/ablab/spades
[sam]: https://en.wikipedia.org/wiki/SAM_(file_format)
[gfa1]: https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
