import { resolve } from 'node:path';
import { defineConfig } from 'vite';

const here = import.meta.dirname;

export default defineConfig({
  // Relative asset URLs, so the built folder can be served from any path.
  base: './',
  // Shared constants live beside the lanes in ../shared.
  server: { fs: { allow: [here, resolve(here, '../shared')] } },
  build: {
    outDir: 'dist',
    emptyOutDir: true,
    chunkSizeWarningLimit: 4000,
    rollupOptions: {
      input: {
        index: resolve(here, 'index.html'),
        shoreline: resolve(here, 'experiences/shoreline/index.html'),
      },
    },
  },
});
