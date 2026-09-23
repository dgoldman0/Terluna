import { resolve } from 'node:path';
import { defineConfig } from 'vite';

const here = import.meta.dirname;

export default defineConfig({
  // Relative asset URLs, so the built folder can be served from any path.
  base: './',
  // The cloud renderer still imports the atmosphere domain's column model from
  // ../atmosphere/column until baked column products replace it.
  server: { fs: { allow: ['..'] } },
  build: {
    outDir: 'dist',
    emptyOutDir: true,
    chunkSizeWarningLimit: 4000,
    rollupOptions: {
      input: {
        index: resolve(here, 'index.html'),
        shoreline: resolve(here, 'experiences/shoreline/index.html'),
        'renderer-lab': resolve(here, 'experiences/renderer-lab/index.html'),
      },
    },
  },
});
