/* Shared registry of engine systems.
 * The shoreline modules were written as browser scripts that attached their
 * systems to one namespace object; modules now import this object explicitly
 * instead of reading a global. The web engine is frozen (see the immersion
 * README), so the registry stays as it is.
 */
export const OM = {};
