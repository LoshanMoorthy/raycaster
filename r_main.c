#define _CRT_SECURE_NO_WARNINGS
#include <SDL.h>
#include <SDL_image.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>

#define CIMGUI_DEFINE_ENUMS_AND_STRUCTS
#include <cimgui.h>

#include "a_assert.h"
#include "r_rand.h"

#define PI 3.14159265359
#define TAU (2.0 * PI)
#define PI_2 (PI / 2.0)
#define PI_4 (PI / 4.0)

#define DEG2RAD(_d) ((_d) * (PI / 180.0))
#define RAD2DEG(_d) ((_d) * (180 / PI))

#define SCREEN_WIDTH 384
#define SCREEN_HEIGHT 216

#define EYE_Z 2
#define HFOV DEG2RAD(90.0)
#define VFOV 0.3

#define ZNEAR 0.0001
#define ZFAR  20.0

typedef struct v2_s { f64 x, y; } v2;
typedef struct v2i_s { i32 x, y; } v2i;

static inline v2i AS_V2I(v2 _v) {
    return (v2i){ .x = _v.x, .y = _v.y };
}

static inline v2 AS_V2(v2i _v) {
    return (v2){ .x = _v.x, .y = _v.y };
}


static inline double dot(v2 _v0, v2 _v1) {
    return (_v0.x * _v1.x) + (_v0.y * _v1.y);
}

static inline double length(v2 _v1) {
    return sqrt(dot(_v1, _v1));
}

static inline v2 normalize(v2 _vn) {
    double l = length(_vn);
    return (v2){ _vn.x / l, _vn.y / l };
}

static inline double cross(v2 _vc0, v2 _vc1) {
    return (_vc0.x * _vc1.y) - (_vc1.x * _vc0.y);
}

// 2D vector rotation
static inline v2 rotate(v2 _vr, double _th) {
    return (v2) {
        (_vr.x * cos(_th)) - (_vr.y * sin(_th)),
        (_vr.x * sin(_th)) + (_vr.y * cos(_th)),
    };
}

// see: https://en.wikipedia.org/wiki/Line–line_intersection

// intersect two infinite lines
static inline v2 intersect_lines(v2 _a0, v2 _a1, v2 _b0, v2 _b1) {
    v2 _l00 = _a0, _l01 = _a1, _l10 = _b0, _l11 = _b1;
    double _d = ((_l00.x - _l01.x) * (_l10.y - _l11.y)) - ((_l00.y - _l01.y) * (_l10.x - _l11.x));
    double _xl0 = cross(_l00, _l01);
    double _xl1 = cross(_l10, _l11);

    return (v2) {
        ((_xl0 * (_l10.x - _l11.x)) - ((_l00.x - _l01.x) * _xl1)) / _d,
            ((_xl0 * (_l10.y - _l11.y)) - ((_l00.y - _l01.y) * _xl1)) / _d
    };
}

// intersect two line segments, returns (nan, nan) if no intersection exists
static inline v2 intersect_segs(v2 _a0, v2 _a1, v2 _b0, v2 _b1) {
    v2 _l00 = _a0, _l01 = _a1, _l10 = _b0, _l11 = _b1;
    double _d = ((_l00.x - _l01.x) * (_l10.y - _l11.y)) - ((_l00.y - _l01.y) * (_l10.x - _l11.x));
    double _t = (((_l00.x - _l10.x) * (_l10.y - _l11.y)) - ((_l00.y - _l10.y) * (_l10.x - _l11.x))) / _d;
    double _u = (((_l00.x - _l10.x) * (_l00.y - _l01.y)) - ((_l00.y - _l10.y) * (_l00.x - _l01.x))) / _d;

    if (fabs(_d) < 0.000001) {
        return (v2) { NAN, NAN };
    }
    else {
        return (_t >= 0 && _t <= 1 && _u >= 0 && _u <= 1) ?
            (v2) {
            _l00.x + (_t * (_l01.x - _l00.x)), _l00.y + (_t * (_l01.y - _l00.y))
        } :
            (v2) {
            NAN, NAN
        };
    }
}

static inline bool point_triangle(v2 _p, v2 _a, v2 _b, v2 _c) {
    double _d = ((_b.y - _c.y) * (_a.x - _c.x) + (_c.x - _b.x) * (_a.y - _c.y));
    double _x = ((_b.y - _c.y) * (_p.x - _c.x) + (_c.x - _b.x) * (_p.y - _c.y)) / _d;
    double _y = ((_c.y - _a.y) * (_p.x - _c.x) + (_a.x - _c.x) * (_p.y - _c.y)) / _d;
    double _z = 1 - _x - _y;

    return (_x > 0) && (_y > 0) && (_z > 0);
}

// -1 right, 0 on, 1 left
static inline int point_side(v2 _p, v2 _a, v2 _b) {
    return -(((_p.x - _a.x) * (_b.y - _a.y)) - ((_p.y - _a.y) * (_b.x - _a.x)));
}

static inline double fract(double _f) {
    return _f - ((i64)(_f));
}

static inline double sign(double _f) {
    return (_f < 0) ? -1.0 : (_f > 0) ? 1.0 : 0.0;
}

static inline double MIN(double _a, double _b) {
    return (_a < _b) ? _a : _b;
}

static inline double MAX(double _a, double _b) {
    return (_a > _b) ? _a : _b;
}

static inline double clamp(double _x, double _mi, double _ma) {
    return min(max(_x, _mi), _ma);
}

static inline void swap(double* _a, double* _b) {
    double _x = *_a;
    *_a = *_b;
    *_b = _x;
}

static inline double ifnan(double _x, double _alt) {
    return isnan(_x) ? _alt : _x;
}

ALWAYS_INLINE void memset32(void* dst, u32 val, usize n) {
    u16* dst32 = dst;
    while (n--) { *dst32++ = val; }
}

ALWAYS_INLINE void memset16(void* dst, u16 val, usize n) {
    u16* dst16 = dst;
    while (n--) { *dst16++ = val; }
}

struct Wall {
    v2i a, b;
    int portal, material;
};

// sector id for "no sector"
#define SECTOR_NONE 0

struct Sector {
    int id;
    usize first_wall, n_walls;

    struct {
        f64 z;
        int texture;
        v2 scale, offset;
    } floor, ceil;
};

enum TextureId {
    TEXTURE_NONE = 0,
    TEXTURE_TEST0,
    TEXTURE_TEST1,
    TEXTURE_TEST2,
    TEXTURE_TEST3,
    TEXTURE_COUNT
};

struct Texture {
    int id;
    u32* data;
    v2i offset, size;
    usize pitch;
};

#define MATERIAL_NONE 0
struct Material {
    int id, texture_id;
    v2 scale;
};

#define TRIG_TAB_SIZE 2048
static f64 sintab[TRIG_TAB_SIZE], tantab[TRIG_TAB_SIZE];

static struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* texture, * debug;
    u32* pixels, * texture_pixels;
    bool quit;

    struct { struct Sector arr[32]; usize n; } sectors;
    struct { struct Wall arr[128]; usize n; } walls;
    struct { struct Material arr[128]; usize n; } materials;
    struct Texture textures[TEXTURE_COUNT];

    struct {
        v2 pos;
        f64 angle, anglecos, anglesin;
        int sector;
    } camera;

    bool sleepy;
} state;

static void make_trigtabs() {
    for (usize i = 0; i < TRIG_TAB_SIZE; i++) {
        // 0..TAU
        const f64 a = (i * TAU) / (f64)(TRIG_TAB_SIZE);
        sintab[i] = sin(a);
        tantab[i] = tan(a);
    }
}

ALWAYS_INLINE f64 fastsin(f64 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return sintab[(TRIG_TAB_SIZE - ((-i) & (TRIG_TAB_SIZE - 1))) & (TRIG_TAB_SIZE - 1)];
    }
    else {
        return sintab[i & (TRIG_TAB_SIZE - 1)];
    }
}


ALWAYS_INLINE f64 fastcos(f64 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return sintab[((-i) + (TRIG_TAB_SIZE / 4)) & (TRIG_TAB_SIZE - 1)];
    }
    else {
        return sintab[(i + (TRIG_TAB_SIZE / 4)) & (TRIG_TAB_SIZE - 1)];
    }
}

ALWAYS_INLINE f64 fasttan(f64 a) {
    const long i = lround(((a) * (TRIG_TAB_SIZE / 2)) / PI);
    if (i < 0) {
        return tantab[(TRIG_TAB_SIZE - ((-i) & (TRIG_TAB_SIZE - 1))) & (TRIG_TAB_SIZE - 1)];
    }
    else {
        return tantab[i & (TRIG_TAB_SIZE - 1)];
    }
}

// stackoverflow.com/questions/42537957
ALWAYS_INLINE f64 fastatan(f64 x) {
#define FT_A 0.0776509570923569
#define FT_B -0.287434475393028
#define FT_C (PI_4 - FT_A - FT_B)
    const f64 xx = x * x;
    return ((FT_A * xx + FT_B) * xx + FT_C) * x;
#undef FT_A
#undef FT_B
#undef FT_C
}

// convert angle in [-(HFOV / 2)..+(HFOV / 2)] to X coordinate relative to
// screen center
ALWAYS_INLINE int screen_angle_to_unbounded_x(f64 angle) {
    // convert to [-PI/4..+PI/4]
    angle = (((angle + (HFOV / 2.0)) / HFOV) * PI_2) - PI_4;
    return (int)((SCREEN_WIDTH / 2) * (1.0 + fasttan(angle)));
}

ALWAYS_INLINE int screen_angle_to_x(f64 angle) {
    return clamp(screen_angle_to_unbounded_x(angle), 0, SCREEN_WIDTH - 1);
}

// convert screen X to [-(HFOV / 2)..+(HFOV / 2)]
ALWAYS_INLINE f64 x_to_screen_angle(int x) {
    // convert back to tan result
    const f64 a = ((x / (f64)(SCREEN_WIDTH / 2)) - 1.0);

    // arctan [-PI/4..PI/4] -> [-1, 1] -> [-HFOV/2..+HFOV/2]
    return (((fastatan(a) + PI_4) / PI_2) * HFOV) - (HFOV / 2.0);
}

// noramlize angle to +/-PI
ALWAYS_INLINE f64 normalize_angle(f64 a) {
    return a - (TAU * floor((a + PI) / TAU));
}

// camera space -> world space (un-rotate and un-translate)
ALWAYS_INLINE v2 camera_pos_to_world(v2 p) {
    const v2
        pr = {
            p.x * state.camera.anglesin + p.y * state.camera.anglecos,
            p.y * state.camera.anglesin - p.x * state.camera.anglecos,
    };
    return (v2) { pr.x + state.camera.pos.x, pr.y + state.camera.pos.y };
}

// convert a wall position AND a screen "y" (floor or ceiling) to a world point
// based on the floor or ceiling height "z"
ALWAYS_INLINE v2 wall_floor_ceiling_to_camera_space(
    v2 dir_to_wall,
    v2 wall_point,
    int y,
    f64 z) {
    const int
        d = (y - (SCREEN_HEIGHT / 2)),
        e = d == 0 ? 1 : d;

    const f64 py = ((VFOV * SCREEN_HEIGHT) * (z - EYE_Z)) / ((f64)e);
    return (v2) {
        wall_point.x - (((wall_point.y - py) / dir_to_wall.y) * dir_to_wall.x),
            py
    };
}

static void present();

// load sectors from file -> state
static int load_sectors(const char* path) {
    printf("Attempting to load sectors from: %s\n", path);

    state.sectors.n = 1;  // Assuming sector 0 does not exist
    state.materials.n = 1;

    FILE* f = fopen(path, "r");
    if (!f) {
        printf("Failed to open the file: %s\n", path);
        return -1;
    }

    int retval = 0;
    enum { SCAN_SECTOR, SCAN_WALL, SCAN_NONE } ss = SCAN_NONE;

    char line[1024], buf[64];
    while (fgets(line, sizeof(line), f)) {
        const char* p = line;
        while (isspace(*p)) {
            p++;
        }

        if (!*p || *p == '#') {
            continue;
        }
        else if (*p == '[') {
            strncpy(buf, p + 1, sizeof(buf));
            const char* section = strtok(buf, "]");
            if (!section) {
                printf("Failed to parse section header in line: %s\n", line);
                retval = -2;
                goto done;
            }

            if (!strcmp(section, "SECTOR")) {
                ss = SCAN_SECTOR;
            }
            else if (!strcmp(section, "WALL")) {
                ss = SCAN_WALL;
            }
            else {
                printf("Unknown section: %s\n", section);
                retval = -3;
                goto done;
            }
        }
        else {
            switch (ss) {
            case SCAN_WALL: {
                struct Wall* wall = &state.walls.arr[state.walls.n++];
                if (sscanf(
                    p,
                    "%d %d %d %d %d",
                    &wall->a.x,
                    &wall->a.y,
                    &wall->b.x,
                    &wall->b.y,
                    &wall->portal) != 5) {
                    printf("Failed to parse wall data: %s\n", p);
                    retval = -4;
                    goto done;
                }
                printf("Parsed wall: a(%d, %d) b(%d, %d) portal(%d)\n",
                       wall->a.x, wall->a.y, wall->b.x, wall->b.y, wall->portal);
                wall->material = 1;
            }; break;
            case SCAN_SECTOR: {
                struct Sector* sector = &state.sectors.arr[state.sectors.n++];
                if (sscanf(
                    p,
                    "%d %" PRIusize " %" PRIusize " %lf %lf",
                    &sector->id,
                    &sector->first_wall,
                    &sector->n_walls,
                    &sector->floor.z,
                    &sector->ceil.z) != 5) {
                    printf("Failed to parse sector data: %s\n", p);
                    retval = -5;
                    goto done;
                }
                printf("Parsed sector: id(%d) first_wall(%zu) n_walls(%zu) floor(%lf) ceil(%lf)\n",
                       sector->id, sector->first_wall, sector->n_walls, sector->floor.z, sector->ceil.z);

                sector->floor.texture = 2;
                sector->floor.scale = (v2){ 1, 1 };
                sector->floor.offset = (v2){ 0, 0 };

                sector->ceil.texture = 3;
                sector->ceil.scale = (v2){ 1, 1 };
                sector->ceil.offset = (v2){ 0, 0 };
            }; break;
            default:
            retval = -6;
            goto done;
            }
        }
    }

    if (ferror(f)) {
        printf("Error reading the file: %s\n", path);
        retval = -128;
        goto done;
    }

done:
    fclose(f);
    return retval;
}

ALWAYS_INLINE u32 texture_sample(
    const struct Texture* tex,
    v2i p) {
    return
        tex->data[
            ((p.y + tex->offset.y) * (tex->pitch / 4))
                + (tex->offset.x + p.x)];
}

static void verline(int x, int y0, int y1, u32 color) {
    for (int y = y0; y <= y1; y++) {
        state.pixels[y * SCREEN_WIDTH + x] = color;
    }
}

static void texline(
    int x,
    int y0,
    int y1,
    int y_floor,
    int y_ceil,
    f64 u,
    const struct Texture* texture) {
    const int tx = (int)(u * texture->size.x);

    const usize pitch_4 = texture->pitch / 4;
    const f64 yd = y_ceil - y_floor;

    for (int y = y0; y <= y1; y++) {
        const int yy = SCREEN_HEIGHT - y - 1;
        const f64 v =
            clamp(1.0 - ((yy - y_floor) / yd), 0.0, 1.0);

        // TODO: TEXTURE_NONE should use pink texture
        const u32 color =
            texture->data[
                ((((int)(v * (texture->size.y - 1))) + texture->offset.y)
                    * pitch_4)
                    + (tx + texture->offset.x)];
        state.pixels[y * SCREEN_WIDTH + x] = color;
    }
}

static void floor_ceil_line(
    int x,
    int y0,
    int y1,
    v2 d_w,
    v2 pw,
    f64 z,
    const struct Texture* tex,
    v2 offset_tex,
    v2 scale_tex) {
    for (int y = y0; y <= y1; y++) {
        const int yy = SCREEN_HEIGHT - y - 1;

        const v2
            pc =
            wall_floor_ceiling_to_camera_space(
                d_w, pw, yy, z),
            pl = camera_pos_to_world(pc);

        // TODO: remove constant
        const v2i pt = {
            (int)(max(pl.x, 0) * 16),
            (int)(max(pl.y, 0) * 16)
        };

        state.pixels[y * SCREEN_WIDTH + x] =
            texture_sample(
                tex,
                (v2i) {
            ((int)(
                (pt.x + offset_tex.x)
                * scale_tex.x))
                % tex->size.x,
                ((int)(
                (pt.y + offset_tex.y)
                * scale_tex.y))
                % tex->size.y
        });
    }
}

// point is in sector if it is on the left side of all walls
static bool point_in_sector(const struct Sector* sector, v2 p) {
    for (usize i = 0; i < sector->n_walls; i++) {
        const struct Wall* wall = &state.walls.arr[sector->first_wall + i];

        if (point_side(p, AS_V2(wall->a), AS_V2(wall->b)) < 0) {
            return false;
        }
    }

    return true;
}

static void render() {
#define FLIP_CLAMP_Y(_y) (clamp(SCREEN_HEIGHT - (_y) - 1, 0, SCREEN_HEIGHT - 1))

    u16 y_top[SCREEN_WIDTH], y_bottom[SCREEN_WIDTH];
    memset16(y_top, SCREEN_HEIGHT - 1, SCREEN_WIDTH);
    memset16(y_bottom, 0, SCREEN_WIDTH);

    // TODO: instead of checking if a sector has been drawn once, count the
    // number of (visible?) portals to the sector and limit based on that
    // instead as there could be multiple portals to the same sector visible
    // number of times each sector has been drawn
    u16 drawcount[ARRLEN(state.sectors.arr)];
    memset16(drawcount, 0, ARRLEN(state.sectors.arr));

    // calculate edges of near/far planes (looking down +Y axis)
    const v2
        zdl = rotate(((v2) { 0.0, 1.0 }), -(HFOV / 2.0)),
        zdr = rotate(((v2) { 0.0, 1.0 }), +(HFOV / 2.0)),
        znl = (v2){ zdl.x * ZNEAR, zdl.y * ZNEAR },
        znr = (v2){ zdr.x * ZNEAR, zdr.y * ZNEAR },
        zfl = (v2){ zdl.x * ZFAR, zdl.y * ZFAR },
        zfr = (v2){ zdr.x * ZFAR, zdr.y * ZFAR };

    enum { QUEUE_MAX = 64 };
    struct QueueEntry { int id, x0, x1; };

    struct { struct QueueEntry arr[QUEUE_MAX]; usize n; } queue = {
        {{ state.camera.sector, 0, SCREEN_WIDTH - 1 }},
        1
    };

    while (queue.n != 0) {
        // grab tail of queue
        struct QueueEntry entry = queue.arr[--queue.n];

        if (drawcount[entry.id] != 0) {
            continue;
        }

        drawcount[entry.id]++;

        const struct Sector* sector = &state.sectors.arr[entry.id];

        const struct Texture
            * tex_floor = &state.textures[sector->floor.texture],
            * tex_ceil = &state.textures[sector->ceil.texture];

        for (usize i = 0; i < sector->n_walls; i++) {
            const struct Wall* wall =
                &state.walls.arr[sector->first_wall + i];

            const struct Material* material =
                &state.materials.arr[wall->material];

            // translate relative to player and rotate points around player's view
            const v2
                qp0 = { wall->a.x - state.camera.pos.x, wall->a.y - state.camera.pos.y },
                qp1 = { wall->b.x - state.camera.pos.x, wall->b.y - state.camera.pos.y },
                op0 = {
                    qp0.x * state.camera.anglesin - qp0.y * state.camera.anglecos,
                    qp0.x * state.camera.anglecos + qp0.y * state.camera.anglesin,
            },
            op1 = {
                qp1.x * state.camera.anglesin - qp1.y * state.camera.anglecos,
                qp1.x * state.camera.anglecos + qp1.y * state.camera.anglesin,
            };

            // wall clipped pos
            v2 cp0 = op0, cp1 = op1;

            // both are negative -> wall is entirely behind player
            if (cp0.y <= 0 && cp1.y <= 0) {
                continue;
            }

            // angle-clip against view frustum
            f64
                ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2),
                ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);

            // clip against view frustum
            const v2
                il = intersect_segs(cp0, cp1, znl, zfl),
                ir = intersect_segs(cp0, cp1, znr, zfr);

            // recompute angles if points change
            if (!isnan(il.x)) {
                cp0 = il;
                ap0 = normalize_angle(atan2(cp0.y, cp0.x) - PI_2);
            }

            if (!isnan(ir.x)) {
                cp1 = ir;
                ap1 = normalize_angle(atan2(cp1.y, cp1.x) - PI_2);
            }

            if (ap0 > ap1) {
                continue;
            }

            if ((ap0 < -(HFOV / 2) && ap1 < -(HFOV / 2))
                || (ap0 > +(HFOV / 2) && ap1 > +(HFOV / 2))) {
                continue;
            }

            // "true" xs before portal clamping
            const int
                tx0 = screen_angle_to_x(ap0),
                tx1 = screen_angle_to_x(ap1);

            // bounds check against portal window
            if (tx0 > entry.x1) { continue; }
            if (tx1 < entry.x0) { continue; }

            const int
                x0 = clamp(tx0, entry.x0, entry.x1),
                x1 = clamp(tx1, entry.x0, entry.x1);

            const f64
                z_floor = sector->floor.z,
                z_ceil = sector->ceil.z,
                nz_floor =
                wall->portal ?
                state.sectors.arr[wall->portal].floor.z : 0,
                nz_ceil =
                wall->portal ?
                state.sectors.arr[wall->portal].ceil.z : 0;

            const f64
                sy0 = ifnan((VFOV * SCREEN_HEIGHT) / cp0.y, 1e10),
                sy1 = ifnan((VFOV * SCREEN_HEIGHT) / cp1.y, 1e10);

            const int
                yf0 = (SCREEN_HEIGHT / 2) + (int)((z_floor - EYE_Z) * sy0),
                yc0 = (SCREEN_HEIGHT / 2) + (int)((z_ceil - EYE_Z) * sy0),
                yf1 = (SCREEN_HEIGHT / 2) + (int)((z_floor - EYE_Z) * sy1),
                yc1 = (SCREEN_HEIGHT / 2) + (int)((z_ceil - EYE_Z) * sy1),
                nyf0 = (SCREEN_HEIGHT / 2) + (int)((nz_floor - EYE_Z) * sy0),
                nyc0 = (SCREEN_HEIGHT / 2) + (int)((nz_ceil - EYE_Z) * sy0),
                nyf1 = (SCREEN_HEIGHT / 2) + (int)((nz_floor - EYE_Z) * sy1),
                nyc1 = (SCREEN_HEIGHT / 2) + (int)((nz_ceil - EYE_Z) * sy1),
                txd = tx1 - tx0,
                yfd = yf1 - yf0,
                ycd = yc1 - yc0,
                nyfd = nyf1 - nyf0,
                nycd = nyc1 - nyc0;

            // compute u boundaries according to wall cutoff
            const f64
                u0 = ((cp0.x - op0.x) / (op1.x - op0.x)),
                u1 = ((cp1.x - op0.x) / (op1.x - op0.x));

            const struct Texture* texture =
                &state.textures[material->texture_id];

            const f64
                u0_z0 = u0 / cp0.y,
                u1_z1 = u1 / cp1.y,
                iz0 = 1.0 / cp0.y,
                iz1 = 1.0 / cp1.y;

            for (int x = x0; x <= x1; x++) {
                // calculate progress along x-axis via tx{0,1} so that walls
                // which are partially cut off due to portal edges still have
                // proper heights
                const f64 xp = ifnan((x - tx0) / (f64)txd, 0);

                // perspective correct texture mapping
                // see en.wikipedia.org/wiki/Texture_mapping
                const f64 ux =
                    clamp(
                        ((((1.0 - xp) * u0_z0) + (xp * u1_z1))
                    / (((1.0 - xp) * iz0) + (xp * iz1))) - 0.0001,
                        0.0,
                        1.0);

                // get y coordinates for this x
                const int
                    tyf = (int)(xp * yfd) + yf0,
                    tyc = (int)(xp * ycd) + yc0,
                    yf = clamp(tyf, y_bottom[x], y_top[x]),
                    yc = clamp(tyc, y_bottom[x], y_top[x]);

                // compute values needed for ceiling/floor texture, which are
                // based off of the fixed world position
                //
                // compute screen angle (theta) and use it to reconstruct:
                // direction to wall (dtw) by inversing the atan2(cp0/1) done
                // above. this will only get us a direction vector, so to find
                // the true point on the wall we need to intersect the line
                // dtw with the original wall in camera space to get us the
                // original wall point(owp)
                const f64 theta = x_to_screen_angle(x);
                const v2
                    d_w = { fastcos(theta + PI_2), fastsin(theta + PI_2) },
                    pw =
                    intersect_lines(
                        ((v2) { 0, 0 }),
                        d_w,
                        cp0,
                        cp1);

                // floor
                if (yf > y_bottom[x]) {
                    floor_ceil_line(
                        x,
                        FLIP_CLAMP_Y(yf),
                        FLIP_CLAMP_Y(y_bottom[x]),
                        d_w,
                        pw,
                        z_floor,
                        tex_floor,
                        sector->floor.offset,
                        sector->floor.scale);
                }

                // ceiling
                if (yc < y_top[x]) {
                    floor_ceil_line(
                        x,
                        FLIP_CLAMP_Y(y_top[x]),
                        FLIP_CLAMP_Y(yc),
                        d_w,
                        pw,
                        z_ceil,
                        tex_ceil,
                        sector->ceil.offset,
                        sector->ceil.scale);
                }

                // flip y-coordinates (SCREEN_HEIGHT - ...) because coordinates
                // are top-left high
                if (wall->portal) {
                    const int
                        tnyf = (int)(xp * nyfd) + nyf0,
                        tnyc = (int)(xp * nycd) + nyc0,
                        nyf = clamp(tnyf, y_bottom[x], y_top[x]),
                        nyc = clamp(tnyc, y_bottom[x], y_top[x]);

                    texline(
                        x,
                        FLIP_CLAMP_Y(yc),
                        FLIP_CLAMP_Y(nyc),
                        tyc,
                        tnyc,
                        ux,
                        &state.textures[2]);

                    texline(
                        x,
                        FLIP_CLAMP_Y(nyf),
                        FLIP_CLAMP_Y(yf),
                        tnyf,
                        tyf,
                        ux,
                        &state.textures[3]);

                    y_top[x] =
                        clamp(
                            min(min(yc, nyc), y_top[x]),
                            0, SCREEN_HEIGHT - 1);

                    y_bottom[x] =
                        clamp(
                            max(max(yf, nyf), y_bottom[x]),
                            0, SCREEN_HEIGHT - 1);
                }
                else {
                    texline(
                        x,
                        FLIP_CLAMP_Y(yc),
                        FLIP_CLAMP_Y(yf),
                        tyf,
                        tyc,
                        ux,
                        texture);

                    y_top[x] = yc;
                    y_bottom[x] = yf;
                }

                if (state.sleepy) {
                    present();
                    SDL_Delay(10);
                }
            }

            if (wall->portal) {
                ASSERT(queue.n != QUEUE_MAX);
                queue.arr[queue.n++] = (struct QueueEntry){
                    .id = wall->portal,
                    .x0 = x0,
                    .x1 = x1
                };
            }
        }
    }

    state.sleepy = false;

#undef FLIP_CLAMP_Y
}

static void present() {
    void* px;
    int pitch;
    SDL_LockTexture(state.texture, NULL, &px, &pitch);
    {
        for (usize y = 0; y < SCREEN_HEIGHT; y++) {
            memcpy(
                &((u8*)px)[y * pitch],
                &state.pixels[y * SCREEN_WIDTH],
                SCREEN_WIDTH * 4);
        }
    }
    SDL_UnlockTexture(state.texture);

    SDL_SetRenderTarget(state.renderer, NULL);
    SDL_SetRenderDrawColor(state.renderer, 0, 0, 0, 0xFF);
    SDL_SetRenderDrawBlendMode(state.renderer, SDL_BLENDMODE_NONE);

    SDL_RenderClear(state.renderer);
    SDL_RenderCopy(state.renderer, state.texture, NULL, NULL);

    SDL_SetTextureBlendMode(state.debug, SDL_BLENDMODE_BLEND);
    SDL_RenderCopy(state.renderer, state.debug, NULL, &((SDL_Rect) { 0, 0, 512, 512 }));
    SDL_RenderPresent(state.renderer);
}

void process_event(SDL_Event* event) {
    ImGuiIO* io = igGetIO();
    switch (event->type) {
    case SDL_MOUSEBUTTONDOWN:
    case SDL_MOUSEBUTTONUP:
    io->MouseDown[event->button.button - 1] = (event->type == SDL_MOUSEBUTTONDOWN);
    break;
    case SDL_MOUSEMOTION:
    io->MousePos.x = (float)event->motion.x;
    io->MousePos.y = (float)event->motion.y;
    break;
    case SDL_MOUSEWHEEL:
    io->MouseWheel += (float)event->wheel.y;
    io->MouseWheelH += (float)event->wheel.x;
    break;
    case SDL_KEYDOWN:
    case SDL_KEYUP: {
        SDL_Keycode key = event->key.keysym.sym;
        bool is_pressed = (event->type == SDL_KEYDOWN);

        ImGuiKey imgui_key;
        switch (key) {
        case SDLK_TAB: imgui_key = ImGuiKey_Tab; break;
        case SDLK_LEFT: imgui_key = ImGuiKey_LeftArrow; break;
        case SDLK_RIGHT: imgui_key = ImGuiKey_RightArrow; break;
        case SDLK_UP: imgui_key = ImGuiKey_UpArrow; break;
        case SDLK_DOWN: imgui_key = ImGuiKey_DownArrow; break;
        case SDLK_PAGEUP: imgui_key = ImGuiKey_PageUp; break;
        case SDLK_PAGEDOWN: imgui_key = ImGuiKey_PageDown; break;
        case SDLK_HOME: imgui_key = ImGuiKey_Home; break;
        case SDLK_END: imgui_key = ImGuiKey_End; break;
        case SDLK_INSERT: imgui_key = ImGuiKey_Insert; break;
        case SDLK_DELETE: imgui_key = ImGuiKey_Delete; break;
        case SDLK_BACKSPACE: imgui_key = ImGuiKey_Backspace; break;
        case SDLK_SPACE: imgui_key = ImGuiKey_Space; break;
        case SDLK_RETURN: imgui_key = ImGuiKey_Enter; break;
        case SDLK_ESCAPE: imgui_key = ImGuiKey_Escape; break;
        default: imgui_key = ImGuiKey_None; break;
        }

        if (imgui_key != ImGuiKey_None) {
            ImGuiIO_AddKeyEvent(io, imgui_key, is_pressed);
        }

        // Handle modifier keys
        io->KeyCtrl = (event->key.keysym.mod & KMOD_CTRL) != 0;
        io->KeyShift = (event->key.keysym.mod & KMOD_SHIFT) != 0;
        io->KeyAlt = (event->key.keysym.mod & KMOD_ALT) != 0;
        io->KeySuper = (event->key.keysym.mod & KMOD_GUI) != 0;
        break;
    }
    case SDL_TEXTINPUT:
    ImGuiIO_AddInputCharactersUTF8(io, event->text.text);
    break;
    }
}




void render_imgui(SDL_Renderer* renderer) {
    ImDrawData* draw_data = igGetDrawData();
    int fb_width = (int)(draw_data->DisplaySize.x * draw_data->FramebufferScale.x);
    int fb_height = (int)(draw_data->DisplaySize.y * draw_data->FramebufferScale.y);

    if (fb_height == 0 || fb_width == 0) return;

    for (int n = 0; n < draw_data->CmdListsCount; n++) {
        ImDrawList* cmd_list = draw_data->CmdLists.Data[n];
        ImDrawVert* vtx_buffer = cmd_list->VtxBuffer.Data;
        ImDrawIdx* idx_buffer = cmd_list->IdxBuffer.Data;

        for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++) {
            ImDrawCmd* pcmd = &cmd_list->CmdBuffer.Data[cmd_i];
            SDL_Rect clip_rect = {
                (int)pcmd->ClipRect.x,
                (int)pcmd->ClipRect.y,
                (int)(pcmd->ClipRect.z - pcmd->ClipRect.x),
                (int)(pcmd->ClipRect.w - pcmd->ClipRect.y)
            };

            SDL_RenderSetClipRect(renderer, &clip_rect);

            for (unsigned int i = 0; i < pcmd->ElemCount; i += 3) {
                SDL_Vertex verts[3];
                for (int j = 0; j < 3; j++) {
                    verts[j].position.x = vtx_buffer[idx_buffer[i + j]].pos.x;
                    verts[j].position.y = vtx_buffer[idx_buffer[i + j]].pos.y;
                    verts[j].color.r = vtx_buffer[idx_buffer[i + j]].col & 0xFF;
                    verts[j].color.g = (vtx_buffer[idx_buffer[i + j]].col >> 8) & 0xFF;
                    verts[j].color.b = (vtx_buffer[idx_buffer[i + j]].col >> 16) & 0xFF;
                    verts[j].color.a = (vtx_buffer[idx_buffer[i + j]].col >> 24) & 0xFF;
                }
                SDL_RenderGeometry(renderer, NULL, verts, 3, NULL, 0);
            }
        }
    }

    SDL_RenderSetClipRect(renderer, NULL);
}


int main(int argc, char* argv[]) {
    make_trigtabs();
    printf("Im here\n");
    int sdlInitResult = SDL_Init(SDL_INIT_VIDEO);
    printf("SDL Init result: %d\n", sdlInitResult);
    ASSERT(!sdlInitResult, "SDL failed to initialize: %s", SDL_GetError());

    int imgInitResult = IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG;
    printf("IMG Init result: %d\n", imgInitResult);
    ASSERT(imgInitResult, "failed to load SDL_image: %s", IMG_GetError());

    state.window =
        SDL_CreateWindow(
            "raycast",
            SDL_WINDOWPOS_UNDEFINED,
            SDL_WINDOWPOS_UNDEFINED,
            1280,
            720,
            0);

    printf("Checking if state.window is valid: %p\n", state.window);
    ASSERT(state.window, "failed to create SDL window: %s\n", SDL_GetError());

    state.renderer =
        SDL_CreateRenderer(
            state.window,
            -1,
            SDL_RENDERER_ACCELERATED
            | SDL_RENDERER_PRESENTVSYNC);

    state.texture =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_STREAMING,
            SCREEN_WIDTH,
            SCREEN_HEIGHT);
    state.debug =
        SDL_CreateTexture(
            state.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_TARGET,
            128,
            128);

    state.pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 4);

    SDL_Surface
        * texture_surface_loaded = IMG_Load("C:/test.png"),
        * texture_surface =
        SDL_ConvertSurfaceFormat(
            texture_surface_loaded,
            SDL_PIXELFORMAT_ABGR8888,
            0);
    state.texture_pixels = texture_surface->pixels;

    state.textures[TEXTURE_TEST0] = (struct Texture){
        .id = TEXTURE_TEST0,
        .data = state.texture_pixels,
        .offset = { 0, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };

    state.textures[TEXTURE_TEST1] = (struct Texture){
        .id = TEXTURE_TEST1,
        .data = state.texture_pixels,
        .offset = { 16, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };


    state.textures[TEXTURE_TEST2] = (struct Texture){
        .id = TEXTURE_TEST2,
        .data = state.texture_pixels,
        .offset = { 32, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };


    state.textures[TEXTURE_TEST3] = (struct Texture){
        .id = TEXTURE_TEST3,
        .data = state.texture_pixels,
        .offset = { 48, 0 },
        .size = { 16, 16 },
        .pitch = texture_surface->pitch
    };

    state.camera.pos = (v2){ 3, 3 };
    state.camera.angle = 0.0;
    state.camera.sector = 1;

    int ret = 0;
    ret = load_sectors("C:/level.txt");
    printf("Load sectors return value: %d\n", ret);
    LOG(
        "loaded %" PRIusize " sectors with %" PRIusize " walls",
        state.sectors.n,
        state.walls.n);
    printf("Checking if sectors were loaded correctly, sectors.n: %zu\n", state.sectors.n);
    ASSERT(state.sectors.n > 1, "No sectors loaded\n");
    // After loading, verify the number of sectors and walls
    state.materials.arr[state.materials.n++] = (struct Material){
        .id = state.materials.n - 1,
        .texture_id = TEXTURE_TEST0,
        .scale = { 1, 1 }
    };

    printf("Im here\n");

    igCreateContext(NULL);
    ImGuiIO* io = igGetIO();

    ImFontAtlas* atlas = io->Fonts;
    ImFontAtlas_AddFontDefault(atlas, NULL);

    // Font atlas
    unsigned char* pixels;
    int width, height, bytes_per_pixel;
    ImFontAtlas_GetTexDataAsRGBA32(atlas, &pixels, &width, &height, &bytes_per_pixel);

    SDL_Texture* font_texture = SDL_CreateTexture(
        state.renderer,
        SDL_PIXELFORMAT_ABGR8888,
        SDL_TEXTUREACCESS_STATIC,
        width,
        height
    );

    SDL_UpdateTexture(font_texture, NULL, pixels, width * 4);
    SDL_SetTextureBlendMode(font_texture, SDL_BLENDMODE_BLEND);

    io->Fonts->TexID = (void*)font_texture;

    while (!state.quit) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            process_event(&ev);
            switch (ev.type) {
            case SDL_QUIT:
            state.quit = true;
            break;
            default:
            break;
            }
        }

        if (state.quit) {
            break;
        }

        const f64 rot_speed = 3.0 * 0.016, move_speed = 3.0 * 0.016;

        const u8* keystate = SDL_GetKeyboardState(NULL);

        if (keystate[SDLK_RIGHT & 0xFFFF]) {
            state.camera.angle += rot_speed;
        }

        if (keystate[SDLK_LEFT & 0xFFFF]) {
            state.camera.angle -= rot_speed;
        }

        state.camera.anglecos = fastcos(state.camera.angle);
        state.camera.anglesin = fastsin(state.camera.angle);

        if (keystate[SDLK_UP & 0xFFFF]) {
            state.camera.pos = (v2){
                state.camera.pos.x + (move_speed * state.camera.anglecos),
                state.camera.pos.y + (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_DOWN & 0xFFFF]) {
            state.camera.pos = (v2){
                state.camera.pos.x - (move_speed * state.camera.anglecos),
                state.camera.pos.y - (move_speed * state.camera.anglesin),
            };
        }

        if (keystate[SDLK_F1 & 0xFFFF]) {
            state.sleepy = true;
        }

        // update player sector
        {
            // BFS neighbors in a circular queue, player is likely to be in one
            // of the neighboring sectors
            enum { QUEUE_MAX = 64 };
            int
                queue[QUEUE_MAX] = { state.camera.sector },
                i = 0,
                n = 1,
                found = SECTOR_NONE;

            while (n != 0) {
                // get front of queue and advance to next
                const int id = queue[i];
                i = (i + 1) % (QUEUE_MAX);
                n--;

                const struct Sector* sector = &state.sectors.arr[id];

                if (point_in_sector(sector, state.camera.pos)) {
                    found = id;
                    break;
                }

                // check neighbors
                for (usize j = 0; j < sector->n_walls; j++) {
                    const struct Wall* wall =
                        &state.walls.arr[sector->first_wall + j];

                    if (wall->portal) {
                        if (n == QUEUE_MAX) {
                            WARN("out of queue space!");
                            goto done;
                        }

                        queue[(i + n) % QUEUE_MAX] = wall->portal;
                        n++;
                    }
                }
            }


        done:
            if (!found) {
                WARN("player is not in a sector!");
                state.camera.sector = 1;
            }
            else {
                state.camera.sector = found;
            }
        }

        memset(state.pixels, 0, SCREEN_WIDTH * SCREEN_HEIGHT * 4);


        int win_width, win_height;
        SDL_GetWindowSize(state.window, &win_width, &win_height);
        io->DisplaySize = (ImVec2){ (float)win_width, (float)win_height };

        igNewFrame();

        igSetNextWindowPos((ImVec2) { 50, 50 }, ImGuiCond_Always, (ImVec2) { 0, 0 });
        igBegin("DEEEBUG", NULL, 0);
        igText("THIS IS A IMGUI TEST");
        igEnd();

        igRender();
        render();
        render_imgui(state.renderer);

        if (!state.sleepy) { present(); }
    }

    igDestroyContext(NULL);
    SDL_FreeSurface(texture_surface_loaded);
    SDL_FreeSurface(texture_surface);
    SDL_DestroyTexture(state.debug);
    SDL_DestroyTexture(state.texture);
    SDL_DestroyRenderer(state.renderer);
    SDL_DestroyWindow(state.window);
    return 0;
}