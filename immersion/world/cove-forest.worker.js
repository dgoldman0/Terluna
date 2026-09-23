/* Places the development cove's far-field trees off the main thread: the same
 * ForestField rule the engine would run. It works through the latest nearest-first
 * queue from the forest, posting each tile as it is made. */
import { OM } from '../engine/om.js';
import { ForestField } from '../engine/forest.js';
import cove from './cove.js';

OM.world = cove;
cove.landscape.initialize();
const field = new ForestField(cove.vegetation, (x, z) => cove.landscape.height(x, z));
let queue = [],
  running = false;
function work() {
  const start = performance.now();
  while (queue.length && performance.now() - start < 25) {
    const [ti, tj] = queue.shift(),
      tile = field.generate(ti, tj);
    self.postMessage({ ti, tj, tile }, [tile.buffer]);
  }
  // Yield so that a newer queue can replace this one.
  running = queue.length > 0;
  if (running) setTimeout(work, 0);
}
self.onmessage = ({ data }) => {
  queue = data.queue.slice();
  if (!running) {
    running = true;
    work();
  }
};
