# visibility-polygon

A library used to construct a visibility polygon for a set of line segments.

## Demo
Original demo by Byron Knoll: http://www.byronknoll.com/visibility.html

## Performance

The time complexity of this implementation is `O(n log n)` (where `n` is the total number of line segments). This is the optimal time complexity for this problem.

## Example Usage

```ts
import {
  breakIntersections,
  convertToSegments,
  compute,
  computeViewport,
  inPolygon,
  Position,
  Polygon,
} from 'visibility-polygon';

const polygons: Polygon[] = [];
// this is the 'world' polygon, which bounds all the polygons you want to compute againts
polygons.push([
  [-1, -1],
  [501, -1],
  [501, 501],
  [-1, 501],
]);
// define vertexes of your polygons
polygons.push([
  [250, 100],
  [260, 140],
  [240, 140],
]);

const segments = breakIntersections(convertToSegments(polygons));

// define your position in which the visibility should be calculated from
const position: Position = [60, 60];

// check if the position is inside the world polygon
if (inPolygon(position, polygons[0])) {
  // compute the visibility polygon, this can be used to draw a polygon with Canvas or WebGL
  const visibility = compute(position, segments);
}
const viewportVisibility = computeViewport(
  position,
  segments,
  [50, 50],
  [450, 450]
);
```

Detailed API information can be found in the usage of the methods in your code-editor (JSDoc).

## Credits

Original source code by Byron Knoll (@byronknoll) on https://github.com/byronknoll/visibility-polygon-js

This version of the library adds TypeScript support and re-implements it in an ESM module, while supporting CommonJS as well.
