/* Shared registry of engine systems (transitional).
 * The shoreline modules were written as browser scripts that attached their
 * systems to one namespace object; the renderer lab swaps material factories
 * through it. Modules now import this object explicitly instead of reading a
 * global. Replacing it with direct imports plus an explicit renderer interface
 * is the next refactor.
 */
export const OM = {};
