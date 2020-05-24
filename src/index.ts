/**
 * Original source code by Byron Knoll (@byronknoll) on https://github.com/byronknoll/visibility-polygon-js
 * This version of the library adds TypeScript support and re-implements it in an ESM module by Liang Chun Wong (@lwong)
 */

export type Position = [number, number];
export type Segments = Position[][];
export type Points = [number, number, number];
export type Polygon = Position[];

/**
 * Computes a visibility polygon. O(N log N) time complexity (where N is the number of line segments).
 * @param position The location of the observer. If the observer is not completely surrounded by line segments, an outer bounding-box will be automatically created (so that the visibility polygon does not extend to infinity).
 * @param segments A list of line segments. Each line segment should be a list of two points. Each point should be a list of two coordinates. Line segments can not intersect each other. Overlapping vertices are OK, but it is not OK if a vertex is touching the middle of a line segment. Use the "breakIntersections" function to fix intersecting line segments.
 * @returns The visibility polygon (in clockwise vertex order).
 */
export function compute(position: Position, segments: Segments) {
  const bounded: Segments = [];
  let minX = position[0];
  let minY = position[1];
  let maxX = position[0];
  let maxY = position[1];
  for (let i = 0; i < segments.length; ++i) {
    for (let j = 0; j < 2; ++j) {
      minX = Math.min(minX, segments[i][j][0]);
      minY = Math.min(minY, segments[i][j][1]);
      maxX = Math.max(maxX, segments[i][j][0]);
      maxY = Math.max(maxY, segments[i][j][1]);
    }
    bounded.push([
      [segments[i][0][0], segments[i][0][1]],
      [segments[i][1][0], segments[i][1][1]],
    ]);
  }
  --minX;
  --minY;
  ++maxX;
  ++maxY;
  bounded.push([
    [minX, minY],
    [maxX, minY],
  ]);
  bounded.push([
    [maxX, minY],
    [maxX, maxY],
  ]);
  bounded.push([
    [maxX, maxY],
    [minX, maxY],
  ]);
  bounded.push([
    [minX, maxY],
    [minX, minY],
  ]);
  let polygon: Position[] = [];
  let sorted = sortPoints(position, bounded);
  let map = new Array(bounded.length);
  for (let i = 0; i < map.length; ++i) map[i] = -1;
  let heap: number[] = [];
  let start: Position = [position[0] + 1, position[1]];
  for (let i = 0; i < bounded.length; ++i) {
    let a1 = angle(bounded[i][0], position);
    let a2 = angle(bounded[i][1], position);
    let active = false;
    if (a1 > -180 && a1 <= 0 && a2 <= 180 && a2 >= 0 && a2 - a1 > 180)
      active = true;
    if (a2 > -180 && a2 <= 0 && a1 <= 180 && a1 >= 0 && a1 - a2 > 180)
      active = true;
    if (active) {
      insert(i, heap, position, bounded, start, map);
    }
  }
  for (let i = 0; i < sorted.length; ) {
    let extend = false;
    let shorten = false;
    let orig = i;
    let vertex = bounded[sorted[i][0]][sorted[i][1]];
    let old_segment = heap[0];
    do {
      if (map[sorted[i][0]] !== -1) {
        if (sorted[i][0] === old_segment) {
          extend = true;
          vertex = bounded[sorted[i][0]][sorted[i][1]];
        }
        remove(map[sorted[i][0]], heap, position, bounded, vertex, map);
      } else {
        insert(sorted[i][0], heap, position, bounded, vertex, map);
        if (heap[0] !== old_segment) {
          shorten = true;
        }
      }
      ++i;
      if (i === sorted.length) break;
    } while (sorted[i][2] < sorted[orig][2] + epsilon());

    if (extend) {
      polygon.push(vertex);
      let cur = intersectLines(
        bounded[heap[0]][0],
        bounded[heap[0]][1],
        position,
        vertex
      );
      if (!equal(cur, vertex)) polygon.push(cur);
    } else if (shorten) {
      polygon.push(
        intersectLines(
          bounded[old_segment][0],
          bounded[old_segment][1],
          position,
          vertex
        )
      );
      polygon.push(
        intersectLines(
          bounded[heap[0]][0],
          bounded[heap[0]][1],
          position,
          vertex
        )
      );
    }
  }
  return polygon;
}

/**
 * Computes a visibility polygon within the given viewport. This can be faster than the "compute" function if there are many segments outside of the viewport.
 * @param position The location of the observer. Must be within the viewport.
 * @param segments A list of line segments. Line segments can not intersect each other. It is OK if line segments intersect the viewport.
 * @param viewportMinCorner The minimum X and Y coordinates of the viewport.
 * @param viewportMaxCorner The maximum X and Y coordinates of the viewport.
 * @returns The visibility polygon within the viewport (in clockwise vertex order).
 */
export function computeViewport(
  position: Position,
  segments: Segments,
  viewportMinCorner: Position,
  viewportMaxCorner: Position
): Position[] {
  let brokenSegments: Position[][] = [];
  let viewport: Position[] = [
    [viewportMinCorner[0], viewportMinCorner[1]],
    [viewportMaxCorner[0], viewportMinCorner[1]],
    [viewportMaxCorner[0], viewportMaxCorner[1]],
    [viewportMinCorner[0], viewportMaxCorner[1]],
  ];
  for (let i = 0; i < segments.length; ++i) {
    if (
      segments[i][0][0] < viewportMinCorner[0] &&
      segments[i][1][0] < viewportMinCorner[0]
    )
      continue;
    if (
      segments[i][0][1] < viewportMinCorner[1] &&
      segments[i][1][1] < viewportMinCorner[1]
    )
      continue;
    if (
      segments[i][0][0] > viewportMaxCorner[0] &&
      segments[i][1][0] > viewportMaxCorner[0]
    )
      continue;
    if (
      segments[i][0][1] > viewportMaxCorner[1] &&
      segments[i][1][1] > viewportMaxCorner[1]
    )
      continue;
    let intersections = [];
    for (let j = 0; j < viewport.length; ++j) {
      let k = j + 1;
      if (k === viewport.length) k = 0;
      if (
        doLineSegmentsIntersect(
          segments[i][0][0],
          segments[i][0][1],
          segments[i][1][0],
          segments[i][1][1],
          viewport[j][0],
          viewport[j][1],
          viewport[k][0],
          viewport[k][1]
        )
      ) {
        let intersect = intersectLines(
          segments[i][0],
          segments[i][1],
          viewport[j],
          viewport[k]
        );
        if (intersect.length !== 2) continue;
        if (
          equal(intersect, segments[i][0]) ||
          equal(intersect, segments[i][1])
        )
          continue;
        intersections.push(intersect);
      }
    }
    let start: Position = [segments[i][0][0], segments[i][0][1]];
    while (intersections.length > 0) {
      let endIndex = 0;
      let endDis = distance(start, intersections[0]);
      for (let j = 1; j < intersections.length; ++j) {
        let dis = distance(start, intersections[j]);
        if (dis < endDis) {
          endDis = dis;
          endIndex = j;
        }
      }
      brokenSegments.push([
        [start[0], start[1]],
        [intersections[endIndex][0], intersections[endIndex][1]],
      ]);
      start[0] = intersections[endIndex][0];
      start[1] = intersections[endIndex][1];
      intersections.splice(endIndex, 1);
    }
    brokenSegments.push([start, [segments[i][1][0], segments[i][1][1]]]);
  }

  let viewportSegments: Segments = [];
  for (let i = 0; i < brokenSegments.length; ++i) {
    if (
      inViewport(brokenSegments[i][0], viewportMinCorner, viewportMaxCorner) &&
      inViewport(brokenSegments[i][1], viewportMinCorner, viewportMaxCorner)
    ) {
      viewportSegments.push([
        [brokenSegments[i][0][0], brokenSegments[i][0][1]],
        [brokenSegments[i][1][0], brokenSegments[i][1][1]],
      ]);
    }
  }
  let eps = epsilon() * 10;
  viewportSegments.push([
    [viewportMinCorner[0] - eps, viewportMinCorner[1] - eps],
    [viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps],
  ]);
  viewportSegments.push([
    [viewportMaxCorner[0] + eps, viewportMinCorner[1] - eps],
    [viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps],
  ]);
  viewportSegments.push([
    [viewportMaxCorner[0] + eps, viewportMaxCorner[1] + eps],
    [viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps],
  ]);
  viewportSegments.push([
    [viewportMinCorner[0] - eps, viewportMaxCorner[1] + eps],
    [viewportMinCorner[0] - eps, viewportMinCorner[1] - eps],
  ]);
  return compute(position, viewportSegments);
}

export function inViewport(
  position: Position,
  viewportMinCorner: Position,
  viewportMaxCorner: Position
): boolean {
  if (position[0] < viewportMinCorner[0] - epsilon()) return false;
  if (position[1] < viewportMinCorner[1] - epsilon()) return false;
  if (position[0] > viewportMaxCorner[0] + epsilon()) return false;
  if (position[1] > viewportMaxCorner[1] + epsilon()) return false;
  return true;
}

/**
 * Calculates whether a point is within a polygon. O(N) time complexity (where N is the number of points in the polygon).
 * @param position The point to check: a list of two coordinates.
 * @param polygon The polygon to check: a list of points. The polygon can be specified in either clockwise or counterclockwise vertex order.
 * @returns True if "position" is within the polygon.
 */
export function inPolygon(position: Position, polygon: Position[]): boolean {
  let val = polygon[0][0];
  for (let i = 0; i < polygon.length; ++i) {
    val = Math.min(polygon[i][0], val);
    val = Math.min(polygon[i][1], val);
  }
  let edge: Position = [val - 1, val - 1];
  let parity = 0;
  for (let i = 0; i < polygon.length; ++i) {
    let j = i + 1;
    if (j === polygon.length) j = 0;
    if (
      doLineSegmentsIntersect(
        edge[0],
        edge[1],
        position[0],
        position[1],
        polygon[i][0],
        polygon[i][1],
        polygon[j][0],
        polygon[j][1]
      )
    ) {
      let intersect = intersectLines(edge, position, polygon[i], polygon[j]);
      if (equal(position, intersect)) return true;
      if (equal(intersect, polygon[i])) {
        if (angle2(position, edge, polygon[j]) < 180) ++parity;
      } else if (equal(intersect, polygon[j])) {
        if (angle2(position, edge, polygon[i]) < 180) ++parity;
      } else {
        ++parity;
      }
    }
  }
  return parity % 2 !== 0;
}

/**
 * Converts the given polygons to list of line segments. O(N) time complexity (where N is the number of polygons).
 * @param polygons a list of polygons (in either clockwise or counterclockwise vertex order). Each polygon should be a list of points. Each point should be a list of two coordinates.
 * @returns a list of line segments.
 */
export function convertToSegments(polygons: Polygon[]): Segments {
  let segments: Segments = [];
  for (let i = 0; i < polygons.length; ++i) {
    for (let j = 0; j < polygons[i].length; ++j) {
      let k = j + 1;
      if (k === polygons[i].length) k = 0;
      segments.push([
        [polygons[i][j][0], polygons[i][j][1]],
        [polygons[i][k][0], polygons[i][k][1]],
      ]);
    }
  }
  return segments;
}

/**
 * Breaks apart line segments so that none of them intersect. O(N^2) time complexity (where N is the number of line segments).
 * @param segments a list of line segments. Each line segment should be a list of two points. Each point should be a list of two coordinates.
 * @returns a list of line segments.
 */
export function breakIntersections(segments: Segments): Segments {
  let output: Segments = [];
  for (let i = 0; i < segments.length; ++i) {
    let intersections = [];
    for (let j = 0; j < segments.length; ++j) {
      if (i === j) continue;
      if (
        doLineSegmentsIntersect(
          segments[i][0][0],
          segments[i][0][1],
          segments[i][1][0],
          segments[i][1][1],
          segments[j][0][0],
          segments[j][0][1],
          segments[j][1][0],
          segments[j][1][1]
        )
      ) {
        let intersect = intersectLines(
          segments[i][0],
          segments[i][1],
          segments[j][0],
          segments[j][1]
        );
        if (intersect.length !== 2) continue;
        if (
          equal(intersect, segments[i][0]) ||
          equal(intersect, segments[i][1])
        )
          continue;
        intersections.push(intersect);
      }
    }
    let start: Position = [segments[i][0][0], segments[i][0][1]];
    while (intersections.length > 0) {
      let endIndex = 0;
      let endDis = distance(start, intersections[0]);
      for (let j = 1; j < intersections.length; ++j) {
        let dis = distance(start, intersections[j]);
        if (dis < endDis) {
          endDis = dis;
          endIndex = j;
        }
      }
      output.push([
        [start[0], start[1]],
        [intersections[endIndex][0], intersections[endIndex][1]],
      ]);
      start[0] = intersections[endIndex][0];
      start[1] = intersections[endIndex][1];
      intersections.splice(endIndex, 1);
    }
    output.push([start, [segments[i][1][0], segments[i][1][1]]]);
  }
  return output;
}

function epsilon(): number {
  return 0.0000001;
}

function equal(a: Position | number[], b: Position | number[]): boolean {
  if (Math.abs(a[0] - b[0]) < epsilon() && Math.abs(a[1] - b[1]) < epsilon())
    return true;
  return false;
}

function remove(
  index: number,
  heap: number[],
  position: Position,
  segments: Segments,
  destination: Position,
  map: Record<number, number>
): void {
  map[heap[index]] = -1;
  if (index === heap.length - 1) {
    heap.pop();
    return;
  }
  heap[index] = heap.pop()!;
  map[heap[index]] = index;
  let cur = index;
  let parent1 = parent(cur);
  if (
    cur !== 0 &&
    lessThan(heap[cur], heap[parent1], position, segments, destination)
  ) {
    while (cur > 0) {
      let parent2 = parent(cur);
      if (
        !lessThan(heap[cur], heap[parent2], position, segments, destination)
      ) {
        break;
      }
      map[heap[parent2]] = cur;
      map[heap[cur]] = parent2;
      let temp = heap[cur];
      heap[cur] = heap[parent2];
      heap[parent2] = temp;
      cur = parent2;
    }
  } else {
    while (true) {
      let left = child(cur);
      let right = left + 1;
      if (
        left < heap.length &&
        lessThan(heap[left], heap[cur], position, segments, destination) &&
        (right === heap.length ||
          lessThan(heap[left], heap[right], position, segments, destination))
      ) {
        map[heap[left]] = cur;
        map[heap[cur]] = left;
        let temp = heap[left];
        heap[left] = heap[cur];
        heap[cur] = temp;
        cur = left;
      } else if (
        right < heap.length &&
        lessThan(heap[right], heap[cur], position, segments, destination)
      ) {
        map[heap[right]] = cur;
        map[heap[cur]] = right;
        let temp = heap[right];
        heap[right] = heap[cur];
        heap[cur] = temp;
        cur = right;
      } else break;
    }
  }
}

function insert(
  index: number,
  heap: number[],
  position: Position,
  segments: Segments,
  destination: Position,
  map: Record<number, number>
): void {
  let intersect = intersectLines(
    segments[index][0],
    segments[index][1],
    position,
    destination
  );
  if ((intersect as Polygon | number[]).length === 0) return;
  let cur = heap.length;
  heap.push(index);
  map[index] = cur;
  while (cur > 0) {
    let parent1 = parent(cur);
    if (!lessThan(heap[cur], heap[parent1], position, segments, destination)) {
      break;
    }
    map[heap[parent1]] = cur;
    map[heap[cur]] = parent1;
    let temp = heap[cur];
    heap[cur] = heap[parent1];
    heap[parent1] = temp;
    cur = parent1;
  }
}

function lessThan(
  index1: number,
  index2: number,
  position: Position,
  segments: Segments,
  destination: Position
): boolean {
  let inter1 = intersectLines(
    segments[index1][0],
    segments[index1][1],
    position,
    destination
  );
  let inter2 = intersectLines(
    segments[index2][0],
    segments[index2][1],
    position,
    destination
  );
  if (!equal(inter1, inter2)) {
    let d1 = distance(inter1 as Position, position);
    let d2 = distance(inter2 as Position, position);
    return d1 < d2;
  }
  let end1 = 0;
  if (equal(inter1, segments[index1][0])) end1 = 1;
  let end2 = 0;
  if (equal(inter2, segments[index2][0])) end2 = 1;
  let a1 = angle2(segments[index1][end1], inter1 as Position, position);
  let a2 = angle2(segments[index2][end2], inter2 as Position, position);
  if (a1 < 180) {
    if (a2 > 180) return true;
    return a2 < a1;
  }
  return a1 < a2;
}

function parent(index: number): number {
  return Math.floor((index - 1) / 2);
}

function child(index: number): number {
  return 2 * index + 1;
}

function angle2(a: Position, b: Position, c: Position): number {
  let a1 = angle(a, b);
  let a2 = angle(b, c);
  let a3 = a1 - a2;
  if (a3 < 0) a3 += 360;
  if (a3 > 360) a3 -= 360;
  return a3;
}

function sortPoints(position: Position, segments: Segments): Points[] {
  let points: Points[] = new Array(segments.length * 2);
  for (let i = 0; i < segments.length; ++i) {
    for (let j = 0; j < 2; ++j) {
      let a = angle(segments[i][j], position);
      points[2 * i + j] = [i, j, a];
    }
  }
  points.sort(function(a, b) {
    return a[2] - b[2];
  });
  return points;
}

function angle(a: Position, b: Position): number {
  return (Math.atan2(b[1] - a[1], b[0] - a[0]) * 180) / Math.PI;
}

function intersectLines(
  a1: Position,
  a2: Position,
  b1: Position,
  b2: Position
): Position {
  let dbx = b2[0] - b1[0];
  let dby = b2[1] - b1[1];
  let dax = a2[0] - a1[0];
  let day = a2[1] - a1[1];

  let u_b = dby * dax - dbx * day;
  if (u_b !== 0) {
    let ua = (dbx * (a1[1] - b1[1]) - dby * (a1[0] - b1[0])) / u_b;
    return [a1[0] - ua * -dax, a1[1] - ua * -day];
  }
  return [] as any;
}

function distance(a: Position | number[], b: Position | number[]): number {
  let dx = a[0] - b[0];
  let dy = a[1] - b[1];
  return dx * dx + dy * dy;
}

function isOnSegment(
  xi: number,
  yi: number,
  xj: number,
  yj: number,
  xk: number,
  yk: number
): boolean {
  return (
    (xi <= xk || xj <= xk) &&
    (xk <= xi || xk <= xj) &&
    (yi <= yk || yj <= yk) &&
    (yk <= yi || yk <= yj)
  );
}

function computeDirection(
  xi: number,
  yi: number,
  xj: number,
  yj: number,
  xk: number,
  yk: number
): -1 | 0 | 1 {
  const a = (xk - xi) * (yj - yi);
  const b = (xj - xi) * (yk - yi);
  return a < b ? -1 : a > b ? 1 : 0;
}

function doLineSegmentsIntersect(
  x1: number,
  y1: number,
  x2: number,
  y2: number,
  x3: number,
  y3: number,
  x4: number,
  y4: number
): boolean {
  const d1 = computeDirection(x3, y3, x4, y4, x1, y1);
  const d2 = computeDirection(x3, y3, x4, y4, x2, y2);
  const d3 = computeDirection(x1, y1, x2, y2, x3, y3);
  const d4 = computeDirection(x1, y1, x2, y2, x4, y4);
  return (
    (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
      ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) ||
    (d1 === 0 && isOnSegment(x3, y3, x4, y4, x1, y1)) ||
    (d2 === 0 && isOnSegment(x3, y3, x4, y4, x2, y2)) ||
    (d3 === 0 && isOnSegment(x1, y1, x2, y2, x3, y3)) ||
    (d4 === 0 && isOnSegment(x1, y1, x2, y2, x4, y4))
  );
}
