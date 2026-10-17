#!/usr/bin/env python3
"""Render the Jura Island flyover: a common-water-level view of LOLA terrain at 16 px/deg, with illustrative lighting.

    python3 visualization/jura-island/render_flyover.py --frames 0 12 22 --scale 0.5   # stills at half size
    python3 visualization/jura-island/render_flyover.py --preview                       # eight stills and a sheet
    python3 visualization/jura-island/render_flyover.py                                 # the 44 s film

Writes to visualization/jura-island/out/ (kept out of Git), or to --out, with render-record.json naming every
product, input and code file it read by sha256. Needs moderngl with an EGL driver; the film also needs ffmpeg.

The land stands at true relief on the reference sphere and the water is the sphere itself, the common level that
every basin below it fills, coloured by depth with 250 m contours. The light comes from a fixed direction in the
north-west, a cartographer's convention, and the haze is a display formula; the illumination models compute neither.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import moderngl
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter

import film
import terrain as terrain_data
from shared.provenance import constants_used
from water_depth import CONTOUR_INTERVAL_M, DEPTH_LEVELS_M, DEPTH_RGB, colours

W, H = film.WIDTH, film.HEIGHT
OUT = HERE / 'out'
SKY = (0.55, 0.68, 0.74)
MAP = dict(south=34.0, north=55.5, west=-51.0, east=-20.0)     # the route map's window (degrees)
MAP_BOX = (W - 284, H - 280, W - 32, H - 62)                   # the route map's place in the frame (pixels)
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
PREVIEW_TIMES = (0, 6, 12, 18, 24, 30, 36, 43.9)
CODE = ('film.py', 'terrain.py', 'water_depth.py', 'render_flyover.py')
VIEW = ('A common-water-level view: LOLA terrain at 16 px/deg less the degree-200 GRAIL geoid, with every basin '
        'below the atlas\'s sea level filled to it, drawn at true relief on the reference sphere with the water '
        'surface on the sphere itself. Depth bands and 250 m contours overlay the water. The light comes from a fixed '
        'direction in the north-west and the haze is a display formula; the illumination models compute neither. '
        'Lake levels, tides, erosion and crustal loading are left out, and travel time is compressed.')

VERTEX = '''
#version 330
in vec3 in_pos;
in vec3 in_normal;
in float in_h;
in vec2 in_uv;
uniform mat4 mvp;
out vec3 pos;
out vec3 normal;
out float height;
out vec2 uv;
void main() {
  pos = in_pos; normal = in_normal; height = in_h; uv = in_uv;
  gl_Position = mvp * vec4(in_pos, 1.0);
}
'''
FRAGMENT = '''
#version 330
in vec3 pos;
in vec3 normal;
in float height;
in vec2 uv;
uniform vec3 eye;
uniform int water;
uniform sampler2D depth_map;
uniform vec3 depth_palette[8];
uniform float depth_edges[7];
uniform float contour_km;
out vec4 frag;
void main() {
  vec3 n = normalize(normal);
  vec3 sun = normalize(vec3(-0.55, 0.65, 0.32));    // illustrative: from the north-west, about 21 degrees up
  vec3 color;
  if (water == 1) {
    float h = height;
    if (uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0) h = texture(depth_map, uv).r;
    float depth = max(-h, 0.0);
    float edge = max(fwidth(depth) * 0.65, 0.0001);
    color = depth_palette[0];
    for (int i = 1; i < 8; i++)
      color = mix(color, depth_palette[i], smoothstep(depth_edges[i - 1] - edge, depth_edges[i - 1] + edge, depth));
    // depth contours, an overlay on the water surface
    float contour_distance = abs(fract(depth / contour_km + 0.5) - 0.5) * contour_km;
    float thickness = max(fwidth(depth) * 0.75, 0.00035);
    float contour = 1.0 - smoothstep(thickness, 2.0 * thickness, contour_distance);
    color = mix(color, vec3(0.82, 0.94, 0.92), 0.38 * contour);
    color *= 0.96 + 0.06 * max(dot(n, sun), 0.0);
  } else {
    float h = clamp(height / 2.8, 0.0, 1.0);
    color = mix(vec3(0.73, 0.67, 0.52), vec3(0.52, 0.51, 0.46), smoothstep(0.0, 0.6, h));
    color = mix(color, vec3(0.85, 0.82, 0.73), smoothstep(0.5, 1.0, h));
    color *= 0.40 + 1.15 * max(dot(n, sun), 0.0);
  }
  float haze = 1.0 - exp(-length(eye - pos) / 1700.0);    // illustrative haze
  color = mix(color, vec3(0.56, 0.68, 0.71), haze * 0.50);
  frag = vec4(pow(clamp(color, 0.0, 1.0), vec3(0.92)), 1.0);
}
'''


def font(size, bold=False):
    return ImageFont.truetype(BOLD if bold else FONT, size)


class Flyover:
    """The scene and its overlays. Overlays are laid out at 1600 x 900 and resized with the frame."""

    def __init__(self, terrain, scale=1.0):
        self.terrain, self.scale = terrain, scale
        self.size = (2 * round(W * scale / 2), 2 * round(H * scale / 2))    # even, as the film's encoder needs
        self.references = terrain_data.references(terrain)
        self.ctx = moderngl.create_standalone_context(backend='egl')
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.disable(moderngl.CULL_FACE)
        self.fbo = self.ctx.simple_framebuffer(self.size, components=3, samples=4)
        self.resolved = self.ctx.simple_framebuffer(self.size, components=3, samples=0)
        self.program = self.ctx.program(vertex_shader=VERTEX, fragment_shader=FRAGMENT)
        self.program['depth_palette'].write(DEPTH_RGB.astype('f4').tobytes())
        self.program['depth_edges'].write((DEPTH_LEVELS_M[1:8] / 1000).astype('f4').tobytes())
        self.program['contour_km'].value = CONTOUR_INTERVAL_M / 1000
        # The depth texture covers the regional crop, its texel centres on the pixel centres
        t = terrain
        texture = self.ctx.texture((t.lon.size, t.lat.size), 1, (t.z / 1000).astype('f4').tobytes(), dtype='f4')
        texture.filter = (moderngl.LINEAR, moderngl.LINEAR)
        texture.repeat_x = texture.repeat_y = False
        texture.use(0)
        self.program['depth_map'].value = 0
        half = 0.5 / t.pixels_per_degree
        self.bounds = (t.lon.min() - half, t.lon.max() + half, t.lat.min() - half, t.lat.max() + half)
        self.meshes = []
        rows = (t.lat >= film.MESH['south']) & (t.lat <= film.MESH['north'])
        cols = (t.lon >= film.MESH['west']) & (t.lon <= film.MESH['east'])
        lat, lon = t.lat[rows], t.lon[cols]
        # A coarser surround at 2 px/deg hides the regional mesh's edge; its land stops 0.25 deg inside that edge
        lon2 = np.append(t.lon4[::2], t.lon4[0] + 360.0)               # the 0/360 seam closed
        z2 = t.z4[::2, ::2] / 1000
        z2 = np.concatenate([z2, z2[:, :1]], axis=1)
        glon, glat = np.meshgrid(lon2, t.lat4[::2])
        cut = (lat.min() + 0.25, lat.max() - 0.25, lon.min() + 0.25, lon.max() - 0.25)
        self.mesh(glon, glat, z2, water=False, cut=cut)
        self.mesh(glon, glat, z2, water=True)
        mlon, mlat = np.meshgrid(lon, lat)
        z = t.z[np.ix_(rows, cols)] / 1000
        self.mesh(mlon, mlat, z, water=False)
        self.mesh(mlon, mlat, z, water=True)
        captions = film.captions(t.sea_level_m, t.water_share, t.pixels_per_degree)
        self.base = self.fit(self.overlay(captions))
        route = self.route_map()
        self.route_at = (round(MAP_BOX[0] * scale), round(MAP_BOX[1] * scale))
        self.route = route if scale == 1 else route.resize(
            (round(route.width * scale), round(route.height * scale)), Image.Resampling.LANCZOS)
        self.chapters = [self.fit(self.chapter_layer(i)) for i in range(len(film.CHAPTERS))]

    def fit(self, layer):
        return layer if layer.size == self.size else layer.resize(self.size, Image.Resampling.LANCZOS)

    def mesh(self, lon, lat, z, water, cut=None):
        if water:
            p = film.xyz(lon, lat, np.zeros_like(z))
            n = film.norm(p + film.R_KM * np.array([0.0, 0.0, 1.0]))
        else:
            p = film.xyz(lon, lat, z)
            # shading normals from heights smoothed slightly at the grid's spacing; the vertices keep the samples
            dr, dc = np.gradient(film.xyz(lon, lat, gaussian_filter(z, 0.65)), axis=(0, 1))
            n = film.norm(np.cross(dr, dc))
        x0, x1, y0, y1 = self.bounds
        uv = np.stack([(lon - x0) / (x1 - x0), (y1 - lat) / (y1 - y0)], axis=-1)
        vertices = np.concatenate([p, n, z[..., None], uv], axis=2).reshape(-1, 9).astype('f4')
        ny, nx = z.shape
        ij = np.arange(ny * nx, dtype=np.uint32).reshape(ny, nx)
        a, b, c, d = ij[:-1, :-1].ravel(), ij[:-1, 1:].ravel(), ij[1:, :-1].ravel(), ij[1:, 1:].ravel()
        indices = np.column_stack([a, c, b, b, c, d])
        if cut is not None:
            south, north, west, east = cut
            mid_lat, mid_lon = (lat[:-1, :-1] + lat[1:, 1:]) / 2, (lon[:-1, :-1] + lon[1:, 1:]) / 2
            inside = (mid_lat > south) & (mid_lat < north) & (mid_lon > west) & (mid_lon < east)
            indices = indices[~inside.ravel()]
        vbo = self.ctx.buffer(vertices.tobytes())
        ibo = self.ctx.buffer(indices.astype('u4').tobytes())
        vao = self.ctx.vertex_array(self.program, [(vbo, '3f 3f 1f 2f', 'in_pos', 'in_normal', 'in_h', 'in_uv')], ibo)
        self.meshes.append((vao, int(water)))

    def overlay(self, captions):
        overlay = Image.new('RGBA', (W, H), (0, 0, 0, 0))
        d = ImageDraw.Draw(overlay)
        # gradients keep the text legible without boxing in the terrain
        for y in range(190):
            d.line((0, y, W, y), fill=(6, 19, 28, int(190 * (1 - y / 190) ** 1.7)))
        for y in range(H - 100, H):
            d.line((0, y, W, y), fill=(5, 17, 23, int(170 * ((y - H + 100) / 100) ** 1.4)))
        d.text((38, 26), 'THE OPEN MOON', font=font(17, True), fill=(209, 228, 225))
        d.text((36, 50), 'JURA ISLAND', font=font(40, True), fill=(248, 246, 233))
        d.text((W - 35, 31), captions['share'], anchor='ra', font=font(18), fill=(240, 244, 239))
        d.text((W - 35, 59), captions['level'], anchor='ra', font=font(16), fill=(200, 219, 219))
        d.text((W - 35, 88), captions['datum'], anchor='ra', font=font(16), fill=(200, 219, 219))
        d.text((36, H - 56), captions['view'], font=font(17), fill=(224, 234, 228))
        d.text((36, H - 31), captions['terrain'], font=font(15), fill=(198, 216, 215))
        d.rounded_rectangle((36, H - 186, 564, H - 80), radius=9, fill=(9, 30, 43, 225))
        d.text((51, H - 178), 'WATER DEPTH BELOW THE SURFACE (m)', font=font(15, True), fill=(230, 240, 230))
        x0, y0, step = 51, H - 147, 60
        for i, c in enumerate(DEPTH_RGB):
            d.rectangle((x0 + i * step, y0, x0 + (i + 1) * step, y0 + 18), fill=tuple(int(v * 255) for v in c) + (255,))
        for i, value in enumerate(DEPTH_LEVELS_M):
            text = f'{int(value):,}' + ('+' if i == len(DEPTH_LEVELS_M) - 1 else '')
            d.text((x0 + i * step, y0 + 23), text, anchor='ma', font=font(12), fill=(220, 236, 231))
        d.rounded_rectangle((36, H - 336, 564, H - 198), radius=9, fill=(9, 30, 43, 225))
        d.text((51, H - 326), 'REFERENCE POINTS  ·  A / B ON ROUTE MAP', font=font(14, True), fill=(204, 225, 223))
        return overlay

    def chapter_layer(self, index):
        """A chapter's fixed panel: its title, two named reference values, and their A and B markers on the route
        map, which stand in place of labels moving over the terrain."""
        layer = Image.new('RGBA', (W, H), (0, 0, 0, 0))
        d = ImageDraw.Draw(layer)
        _, _, title, names = film.CHAPTERS[index]
        d.text((38, 108), title, font=font(26), fill=(233, 236, 220))
        for row, name in enumerate(names):
            lon, lat, value, sense = self.references[name]
            letter, cy = 'AB'[row], H - 278 + row * 43
            d.ellipse((51, cy - 12, 75, cy + 12), fill=(223, 232, 209))
            d.text((63, cy), letter, font=font(14, True), anchor='mm', fill=(12, 36, 45))
            d.text((88, cy), name, font=font(18, True), anchor='lm', fill=(234, 239, 226))
            d.text((546, cy), f'{value:,} m {sense}', font=font(18, True), anchor='rm', fill=(255, 239, 192))
            x, y = self.map_point(lon, lat)
            shift = -12 if row == 0 else 12
            lx, ly = x + shift, y + shift
            d.line((x, y, lx, ly), fill=(235, 238, 220), width=2)
            d.ellipse((x - 2, y - 2, x + 2, y + 2), fill=(235, 238, 220))
            d.ellipse((lx - 9, ly - 9, lx + 9, ly + 9), fill=(223, 232, 209), outline=(12, 36, 45), width=1)
            d.text((lx, ly), letter, font=font(12, True), anchor='mm', fill=(12, 36, 45))
        return layer

    def route_map(self):
        t = self.terrain
        rows = (t.lat >= MAP['south']) & (t.lat <= MAP['north'])
        cols = (t.lon >= MAP['west']) & (t.lon <= MAP['east'])
        z = t.z[np.ix_(rows, cols)]
        dy, dx = np.gradient(gaussian_filter(z, 1))
        relief = np.clip(0.87 + (0.5 * dx - 0.5 * dy) / 800, 0.63, 1.12)
        colour = np.where((z > 0)[..., None], np.array([0.74, 0.68, 0.53]) * relief[..., None],
                          colours(np.maximum(-z, 0)))
        width, height = MAP_BOX[2] - MAP_BOX[0], MAP_BOX[3] - MAP_BOX[1]
        im = Image.fromarray((np.clip(colour, 0, 1) * 255).astype('uint8')).resize(
            (width, height), Image.Resampling.LANCZOS).convert('RGBA')
        md = ImageDraw.Draw(im)
        md.rectangle((0, 0, width - 1, height - 1), outline=(197, 222, 215, 255), width=1)
        md.text((236, 9), 'N', anchor='ra', font=font(13, True), fill='white')
        md.line((230, 39, 230, 27), fill='white', width=1)
        md.polygon([(230, 22), (226, 29), (234, 29)], fill='white')
        md.text((8, 8), 'ROUTE', font=font(12, True), fill=(229, 236, 222))
        md.line([self.map_point(*film.path(s)[:2], local=True) for s in np.linspace(0, film.DURATION, 180)],
                fill=(228, 231, 205, 170), width=1)
        km = 100 / (film.R_KM * np.pi / 180 * np.cos(np.radians(45))) / (MAP['east'] - MAP['west']) * width
        md.line((12, 199, 12 + km, 199), fill=(234, 239, 221), width=2)
        md.text((12, 180), '100 km', font=font(12), fill=(234, 239, 221))
        return im

    def map_point(self, lon, lat, local=False):
        x = (lon - MAP['west']) / (MAP['east'] - MAP['west']) * (MAP_BOX[2] - MAP_BOX[0])
        y = (MAP['north'] - lat) / (MAP['north'] - MAP['south']) * (MAP_BOX[3] - MAP_BOX[1])
        return (x, y) if local else (MAP_BOX[0] + x, MAP_BOX[1] + y)

    def frame(self, t):
        t = float(t)
        lon, lat, alt = film.path(t)[:3]
        eye, view = film.camera_pose(t)
        projection = film.perspective(film.FIELD_OF_VIEW_DEG, W / H, max(0.5, alt * 0.10), 4000.0)
        self.program['mvp'].write((projection @ view).T.astype('f4').tobytes())
        self.program['eye'].value = tuple(float(v) for v in eye)
        self.fbo.use()
        self.fbo.clear(*SKY, 1.0)
        for vao, water in self.meshes:
            self.program['water'].value = water
            vao.render(moderngl.TRIANGLES)
        self.ctx.copy_framebuffer(self.resolved, self.fbo)
        w, h = self.size
        rgb = np.frombuffer(self.resolved.read(components=3, alignment=1), dtype=np.uint8).reshape(h, w, 3)[::-1]
        im = Image.fromarray(rgb).convert('RGBA')
        im.alpha_composite(self.base)
        im.alpha_composite(self.route, self.route_at)
        chapter, opacity = film.chapter_state(t)
        layer = self.chapters[chapter]
        if opacity < 1:
            layer = layer.copy()
            layer.putalpha(layer.getchannel('A').point(lambda a: round(a * opacity)))
        im.alpha_composite(layer)
        d, s = ImageDraw.Draw(im), self.scale
        # the position marker moves smoothly and has no heading arrow to spin when the camera looks down
        x, y = self.map_point(lon, lat)
        if MAP_BOX[0] < x < MAP_BOX[2] and MAP_BOX[1] < y < MAP_BOX[3]:
            d.ellipse(((x - 5) * s, (y - 5) * s, (x + 5) * s, (y + 5) * s), fill=(255, 217, 123),
                      outline=(16, 36, 40), width=1)
        d.line((36 * s, (H - 5) * s, (36 + (W - 72) * t / film.DURATION) * s, (H - 5) * s), fill=(229, 210, 153),
               width=max(1, round(3 * s)))
        return im.convert('RGB')


def contact_sheet(stills):
    """The stills at 800 x 450 in two columns, each marked with its time."""
    rows = (len(stills) + 1) // 2
    sheet = Image.new('RGB', (1600, rows * 450), (7, 20, 26))
    for i, (t, image) in enumerate(stills):
        thumb = image.resize((800, 450), Image.Resampling.LANCZOS)
        ImageDraw.Draw(thumb).text((10, 420), f'{t:.1f} s', font=font(18), fill='white')
        sheet.paste(thumb, ((i % 2) * 800, (i // 2) * 450))
    return sheet


def encode(movie, path):
    """The whole film through ffmpeg's H.264 encoder."""
    w, h = movie.size
    t = movie.terrain
    captions = film.captions(t.sea_level_m, t.water_share, t.pixels_per_degree)
    command = ['ffmpeg', '-y', '-v', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', f'{w}x{h}',
               '-r', str(film.FPS), '-i', '-', '-an', '-c:v', 'libx264', '-preset', 'medium', '-crf', '19',
               '-pix_fmt', 'yuv420p', '-movflags', '+faststart', '-metadata', 'title=Jura Island flyover',
               '-metadata', f'comment={captions["view"]}; {captions["terrain"]}', str(path)]
    process = subprocess.Popen(command, stdin=subprocess.PIPE)
    start = time.monotonic()
    for i, s in enumerate(film.frame_times()):
        process.stdin.write(movie.frame(s).tobytes())
        if i % (film.FPS * 2) == 0:
            print(f'{s:.0f}/{film.DURATION} s rendered; {time.monotonic() - start:.1f} s elapsed', flush=True)
    process.stdin.close()
    if process.wait() != 0:
        raise RuntimeError('Video encoding failed')


def render_record(data, movie, written):
    files = [HERE / f for f in CODE] + [ROOT / 'geography/topography.py', ROOT / 'geography/fetch_inputs.py']
    area = data.cell_area_km2()
    return dict(schema='terluna.visualization.jura-flyover-render/1', view=VIEW,
                code={terrain_data.relative(p): terrain_data.digest(p) for p in files},
                constants=constants_used(files), read=data.hashes,
                sea_level_m=data.sea_level_m, water_share=data.water_share,
                island=dict(area_km2=round(float(area[data.island].sum()), 2),
                            summit_m=round(float(data.z[data.island].max()), 2)),
                renderer=movie.ctx.info['GL_RENDERER'], size=list(movie.size), frame_rate=film.FPS,
                written=written, camera=film.camera_checks(), soundings=terrain_data.soundings(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--frames', type=float, nargs='+', help='render stills at these times (s)')
    parser.add_argument('--preview', action='store_true', help='render eight stills and a contact sheet')
    parser.add_argument('--scale', type=float, default=1.0, help='frame size as a fraction of 1600 x 900')
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    data = terrain_data.build()
    movie = Flyover(data, args.scale)
    print('Renderer:', movie.ctx.info['GL_RENDERER'], flush=True)
    written = []
    if args.frames or args.preview:
        stills = []
        for t in args.frames or PREVIEW_TIMES:
            image = movie.frame(t)
            name = f'frame-{t:04.1f}s.png'
            image.save(args.out / name)
            written.append(name)
            stills.append((t, image))
        if args.preview:
            contact_sheet(stills).save(args.out / 'contact-sheet.jpg', quality=90)
            written.append('contact-sheet.jpg')
    else:
        encode(movie, args.out / 'jura-island-flyover.mp4')
        written.append('jura-island-flyover.mp4')
    (args.out / 'render-record.json').write_text(json.dumps(render_record(data, movie, written), indent=1) + '\n')
    print(args.out, flush=True)


if __name__ == '__main__':
    main()
