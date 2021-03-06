Release 0.5-r234 (21 March 2021)
--------------------------------

Many changes:

 * Use Tarjan's SCC algorithm to identify bubble. This algorithm is more robust
   to inversions, but requires reference paths.

 * Generate unitig sequences from assembly graphs.

 * Use BFS to extract subgraphs within a radius. The older code uses DFS which
   is sometimes counter-intuitive.

(21 March 2021, r234)



Release 0.4-r165 (23 December 2019)
-----------------------------------

Many changes since the last release:

 * Revamped assembly graph cleaning. Now most cleaning operations are both
   count- and length-aware. Bubble popping become more precise. Also added
   topology-aware edge cutting that cuts weak edges while preserving the
   topological contiguity.

 * Added the paf2gfa tool. It is essentially minisam combined with the new
   gfatools graph cleaning APIs. Currently it is specifically designed for PAFs
   generated by hifiasm.

 * More meaningful graph linearization over inversions. The "blacklist" command
   now ignores inversions.

 * Subgraph extraction now considers both directions to properly deal with
   starting segments in bubbles.

 * Added the "bubble" command to output information of each bubble.

(23 December 2019, r165)



Release 0.3-r70 (22 August 2019)
--------------------------------

Notable changes:

 * Generate blacklisted regions from rGFA. The current implementation assumes a
   minigraph graph. It is not general.

 * Output basic statistics of GFA graphs.

 * Check and remove multiple edges.

(22 August 2019, r70)



Release 0.2-r51 (7 August 2019)
-------------------------------

Support SR tags on L-lines.

(7 August 2019, r51)



Release 0.1-r40 (6 July 2019)
-----------------------------

Initial release.

(6 July 2019, r40)
