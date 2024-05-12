/** @type {HTMLCanvasElement}*/
const canvas = document.querySelector("canvas");
/** @type {CanvasRenderingContext2D} */
const ctx = canvas.getContext("2d");

windowRisized();
window.addEventListener("resize", (e) => {
    windowRisized();
});

function windowRisized() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
}

class Vector3D {
    /**
     *
     * @param {number} x
     * @param {number} y
     * @param {number} z
     */
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {Vector3D}
     */
    static add(v1, v2) {
        return new Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {Vector3D}
     */
    static sub(v1, v2) {
        return new Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {number} s
     * @returns {Vector3D}
     */
    static multiply(v1, s) {
        return new Vector3D(v1.x * s, v1.y * s, v1.z * s);
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {number} s
     * @returns {Vector3D}
     */
    static divide(v1, s) {
        return new Vector3D(v1.x / s, v1.y / s, v1.z / s);
    }
    /**
     *
     * @param {Vector3D} v
     * @returns {Vector3D}
     */
    static magnitude(v) {
        return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }
    /**
     *
     * @param {Vector3D} v
     * @returns {Vector3D}
     */
    static squaredMagnitude(v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }
    /**
     *
     * @param {Vector3D} v
     * @returns {Vector3D}
     */
    static normalize(v) {
        const mag = Vector3D.magnitude(v);
        return new Vector3D(v.x / mag, v.y / mag, v.z / mag);
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {number}
     */
    static dot(v1, v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {Vector3D}
     */
    static cross(v1, v2) {
        return new Vector3D(
            v1.y * v2.z - v1.z * v2.y,
            -1 * (v1.x * v2.z - v1.z * v2.x),
            v1.x * v2.y - v1.y * v2.x
        );
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {Vector3D}
     */
    static lerp(v1, v2, t) {
        return Vector3D.add(v1, Vector3D.multiply(Vector3D.sub(v2, v1), t));
    }
    /**
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @returns {number}
     * */
    static distance(v1, v2) {
        return Vector3D.magnitude(Vector3D.sub(v1, v2));
    }
}
class Plane {
    /**
     *
     * @param {number} a
     * @param {number} b
     * @param {number} c
     * @param {number} d
     */
    constructor(ax, by, cz, d) {
        this.ax = ax;
        this.by = by;
        this.cz = cz;
        this.d = d;
    }
    /**
     *
     * @param {Vector3D} v1
     * @param {Vector3D} v2
     * @param {Vector3D} v3
     * vectors are considered to start from the same point
     * @returns {Plane}
     */
    static planeFromVectors(v1, v2, v3) {
        let a = Vector3D.sub(v2, v1);
        let b = Vector3D.sub(v3, v1);
        let n = Vector3D.normalize(Vector3D.cross(a, b));
        let d = -n.x * v1.x - n.y * v1.y - n.z * v1.z;
        return new Plane(n.x, n.y, n.z, d);
    }
    /**
     *
     * @param {Plane} plane
     * @returns {Vector3D}
     */
    static normal(plane) {
        return new Vector3D(plane.ax, plane.by, plane.cz);
    }
    /**
     *
     * @param {Ray} ray
     * @param {Plane} plane
     * @returns {number} parametric value of the ray at the intersection point
     */
    static planeRayIntersection(ray, plane) {
        let t =
            -(
                plane.ax * ray.start.x +
                plane.by * ray.start.y +
                plane.cz * ray.start.z +
                plane.d
            ) /
            (plane.ax * (ray.end.x - ray.start.x) +
                plane.by * (ray.end.y - ray.start.y) +
                plane.cz * (ray.end.z - ray.start.z));
        return t;
    }
}
class Ray {
    /**
     *
     * @param {Vector3D} start
     * @param {Vector3D} end
     * @param {string} color
     */
    constructor(start, end, color) {
        this.start = start;
        this.end = end;
        this.color = color ? color : "white";
    }
    /**
     *
     * @param {Ray} ray
     * @param {number} t
     * @returns {Vector3D}
     */
    static pointAtParameter(ray, t) {
        return Vector3D.lerp(ray.start, ray.end, t);
    }
    /**
     *
     * @param {Ray} ray
     * @param {Plane} plane
     * @returns {number}
     */
    static planeRayIntersection(ray, plane) {
        let t =
            -(
                plane.ax * ray.start.x +
                plane.by * ray.start.y +
                plane.cz * ray.start.z +
                plane.d
            ) /
            (plane.ax * (ray.end.x - ray.start.x) +
                plane.by * (ray.end.y - ray.start.y) +
                plane.cz * (ray.end.z - ray.start.z));
        return t;
    }
}
class Triangle {
    /**
     *
     * @param {Vector3D} a
     * @param {Vector3D} b
     * @param {Vector3D} c
     * @param {string} color
     */
    constructor(a, b, c, color) {
        this.a = a;
        this.b = b;
        this.c = c;
        this.color = color ? color : "white";
        this.plane = Plane.planeFromVectors(a, b, c);
        let la = Vector3D.distance(a, b);
        let lb = Vector3D.distance(b, c);
        let lc = Vector3D.distance(c, a);
        this.perimeter = la + lb + lc;
        let s = this.perimeter / 2;
        this.area = Math.sqrt(s * (s - la) * (s - lb) * (s - lc));
    }
    /**
     *
     * @param {Ray} ray
     * @param {Triangle} triangle
     * @returns {boolean | number} false if no intersection,
     * parametric value of the ray at the intersection point otherwise
     */
    static triangleRayIntersection(ray, triangle) {
        let t = Ray.planeRayIntersection(ray, triangle.plane);
        let i = Ray.pointAtParameter(ray, t);
        let t1 = new Triangle(i, triangle.a, triangle.b);
        let t2 = new Triangle(i, triangle.b, triangle.c);
        let t3 = new Triangle(i, triangle.c, triangle.a);
        if (
            (t1.area + t2.area + t3.area) * 1000000000 <=
            triangle.area + triangle.area * 1000000000
        ) {
            return t;
        }
        return false;
    }
}
class MatrixVector4D {
    /**
     *
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @param {number} w is the translation
     */
    constructor(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
    /**
     *
     * @param {number} s
     * @returns {MatrixVector4D}
     */
    scale(s) {
        return new MatrixVector4D(this.x * s, this.y * s, this.z * s, this.w);
    }
}
class Matrix {
    /**
     *
     * @param {MatrixVector4D} i
     * @param {MatrixVector4D} j
     * @param {MatrixVector4D} k
     * by default vectors are the basis of the space with no translation
     */
    constructor(i, j, k) {
        this.i = i ? i : new MatrixVector4D(1, 0, 0, 0);
        this.j = j ? j : new MatrixVector4D(0, 1, 0, 0);
        this.k = k ? k : new MatrixVector4D(0, 0, 1, 0);
    }
    /**
     *
     * @param {Matrix} m
     * @param {number} s
     * @returns {Matrix}
     */
    static scale(m, s) {
        let newM = new Matrix(m.i.scale(s), m.j.scale(s), m.k.scale(s));
        return newM;
    }
    /**
     *
     * @param {Matrix} m
     * @param {Vector3D} v
     * @returns {Vector3D}
     */
    static multiplyVector(m, v) {
        let newV = new Vector3D(
            m.i.x * v.x + m.j.x * v.y + m.k.x * v.z + m.i.w,
            m.i.y * v.x + m.j.y * v.y + m.k.y * v.z + m.j.w,
            m.i.z * v.x + m.j.z * v.y + m.k.z * v.z + m.k.w
        );
        return newV;
    }
    /**
     *
     * @param {Matrix} m2
     * @param {Matrix} m1
     * @returns {Matrix}
     */
    static multiplyMatrix(m1, m2) {
        let i = new MatrixVector4D(
            m1.i.x * m2.i.x + m1.j.x * m2.i.y + m1.k.x * m2.i.z,
            m1.i.y * m2.i.x + m1.j.y * m2.i.y + m1.k.y * m2.i.z,
            m1.i.z * m2.i.x + m1.j.z * m2.i.y + m1.k.z * m2.i.z,
            m1.i.w
        );
        let j = new MatrixVector4D(
            m1.i.x * m2.j.x + m1.j.x * m2.j.y + m1.k.x * m2.j.z,
            m1.i.y * m2.j.x + m1.j.y * m2.j.y + m1.k.y * m2.j.z,
            m1.i.z * m2.j.x + m1.j.z * m2.j.y + m1.k.z * m2.j.z,
            m1.j.w
        );
        let k = new MatrixVector4D(
            m1.i.x * m2.k.x + m1.j.x * m2.k.y + m1.k.x * m2.k.z,
            m1.i.y * m2.k.x + m1.j.y * m2.k.y + m1.k.y * m2.k.z,
            m1.i.z * m2.k.x + m1.j.z * m2.k.y + m1.k.z * m2.k.z,
            m1.k.w
        );
        return new Matrix(i, j, k);
    }
    /**
     * @param {Matrix} m
     * @param {number} angle
     * @returns {Matrix}
     * rotate the matrix around the x(i) axis
     * rotate the yz plane
     * */
    static rotateX(m, angle) {
        let cos = Math.cos(angle);
        let sin = Math.sin(angle);
        let i = new MatrixVector4D(1, 0, 0, 0);
        let j = new MatrixVector4D(0, cos, -sin, 0);
        let k = new MatrixVector4D(0, sin, cos, 0);
        return Matrix.multiplyMatrix(m, new Matrix(i, j, k));
    }
    /**
     * @param {Matrix} m
     * @param {number} angle
     * @returns {Matrix}
     * rotate the matrix around the y(j) axis
     * rotate the xz plane
     * */
    static rotateY(m, angle) {
        let cos = Math.cos(angle);
        let sin = Math.sin(angle);
        let i = new MatrixVector4D(cos, 0, sin, 0);
        let j = new MatrixVector4D(0, 1, 0, 0);
        let k = new MatrixVector4D(-sin, 0, cos, 0);
        return Matrix.multiplyMatrix(m, new Matrix(i, j, k));
    }
    /**
     * @param {Matrix} m
     * @param {number} angle
     * @returns {Matrix}
     * rotate the matrix around the z(k) axis
     * rotate the xy plane
     * */
    static rotateZ(m, angle) {
        let cos = Math.cos(angle);
        let sin = Math.sin(angle);
        let i = new MatrixVector4D(cos, -sin, 0, 0);
        let j = new MatrixVector4D(sin, cos, 0, 0);
        let k = new MatrixVector4D(0, 0, 1, 0);
        return Matrix.multiplyMatrix(m, new Matrix(i, j, k));
    }
}

class Film {
    /**
     * @param {Vector3D} position is just the traslation of the film
     * @param {number} width
     * @param {number} height
     * @param {CanvasRenderingContext2D} ctx
     * @param {number} pixelSize
     * Film is defined by a position and a width and height
     * recommended position is (0, 0, 0)
     */
    constructor(position, width, height, ctx, pixelSize) {
        this.position = position;
        this.width = width;
        this.height = height;
        this.halfWidth = width / 2;
        this.halfHeight = height / 2;
        this.ctx = ctx;
        this.pixelSize = pixelSize;
        this.pixelsOnFilm = [];
        this.pixelsInSpace = [];
        /**
         * creating the pixels on the film and
         * their corresponding position in the space
         * pixel goes from left to right row by row
         */
        for (
            let i = this.halfHeight;
            i >= -this.halfHeight;
            i = i - pixelSize
        ) {
            for (
                let j = -this.halfWidth;
                j <= this.halfWidth;
                j = j + pixelSize
            ) {
                this.pixelsOnFilm.push({
                    x: j + this.halfWidth,
                    y: -i + this.halfHeight,
                });
                this.pixelsInSpace.push(new Vector3D(j, i, position.z));
            }
        }
    }
    /**
     * x and y are the position of the pixel on the film
     * @param {number} x
     * @param {number} y
     * @param {string} color
     */
    drawPixel(x, y, color) {
        this.ctx.fillStyle = color;
        this.ctx.fillRect(x, y, this.pixelSize, this.pixelSize);
    }
    /**
     * @param {number} i is the index of the pixel on the film
     * @param {string} color
     */
    drawPixelbyIndex(i, color) {
        this.ctx.fillStyle = color;
        this.ctx.fillRect(
            this.pixelsOnFilm[i].x,
            this.pixelsOnFilm[i].y,
            this.pixelSize,
            this.pixelSize
        );
    }
}
class Camera {
    /**
     * @param {Film} film
     * @param {number} distance of film from the camera (fov)
     * @param {Matrix} matrix
     * */
    constructor(film, distance, matrix) {
        this.matrix = matrix;
        if (!matrix) {
            this.matrix = new Matrix(
                new MatrixVector4D(1, 0, 0, 0),
                new MatrixVector4D(0, 1, 0, 0),
                new MatrixVector4D(0, 0, 1, 0)
            );
        }
        this.film = film;
        this.distance = distance;
        /**
         * pixels in the space that the camera will trace through
         */
        this.traceThroughPixels = [];
        this.film.pixelsInSpace.forEach((p, i) => {
            this.traceThroughPixels[i] = Matrix.multiplyVector(
                this.matrix,
                new Vector3D(p.x, p.y, this.distance)
            );
        });
    }
    /**
     *
     * @param {number} i is the index of the pixel on the film
     * @returns {Ray} from the camera to the pixel in the space
     */
    trace(i) {
        let start = new Vector3D(
            this.matrix.i.w,
            this.matrix.j.w,
            this.matrix.k.w
        );
        let end = this.traceThroughPixels[i];
        return new Ray(start, end);
    }
    /**
     *
     * @param {Matrix} matrix
     */
    setMatrix(matrix) {
        this.matrix = matrix;
        this.film.pixelsInSpace.forEach((p, i) => {
            this.traceThroughPixels[i] = Matrix.multiplyVector(
                this.matrix,
                new Vector3D(p.x, p.y, this.distance)
            );
        });
    }
}

/**
 * @type {Array<Mesh>}
 */
const meshes = [];
class Mesh {
    /**
     * @param {string} file
     * @param {Matrix} matrix
     * file should be a json, and structed as follows
     * [{"vertices":[[0,0,0],[0,0,0],[0,0,0]], "color":"rgba(0,0,0,1)"}]
     * if no matrix is provided the mesh will be at 0,0,1000
     */
    constructor(file, matrix) {
        this.matrix = matrix;
        if (!matrix) {
            this.matrix = new Matrix(
                new MatrixVector4D(1, 0, 0, 0),
                new MatrixVector4D(0, 1, 0, 0),
                new MatrixVector4D(0, 0, 1, 1000)
            );
        }
        this.originalTriangles = [];
        this.triangles = [];
        this.load(file)
            .then((triangles) => {
                this.originalTriangles = triangles;
                this.triangles = triangles;
                this.alignTriangles();
                meshes.push(this);
            })
            .catch((e) => {
                console.error(e);
            });
    }
    /**
     * @param {string} file
     */
    async load(file) {
        let response = await fetch(file);
        let data = await response.json();
        let triangles = data.map((t) => {
            const v = t.vertices;
            return new Triangle(
                new Vector3D(v[0][0], v[0][1], v[0][2]),
                new Vector3D(v[1][0], v[1][1], v[1][2]),
                new Vector3D(v[2][0], v[2][1], v[2][2]),
                t.color
            );
        });
        return triangles;
    }
    /**
     * align the triangles with the matrix.
     * to call after any transformation on the mesh
     * by not calling Mesh methods
     */
    alignTriangles() {
        this.triangles = this.originalTriangles.map((t) => {
            return new Triangle(
                Matrix.multiplyVector(this.matrix, t.a),
                Matrix.multiplyVector(this.matrix, t.b),
                Matrix.multiplyVector(this.matrix, t.c),
                t.color
            );
        });
    }
    /**
     * numbers are in radians
     * @param {number} alpha
     * @param {number} beta
     * @param {number} theta
     */
    rotateXYZ(alpha, beta, theta) {
        this.matrix = Matrix.rotateX(this.matrix, alpha);
        this.matrix = Matrix.rotateY(this.matrix, beta);
        this.matrix = Matrix.rotateZ(this.matrix, theta);
        this.alignTriangles();
    }
    /**
     *
     * @param {number} x
     * @param {number} y
     * @param {number} z
     */
    translate(x, y, z) {
        this.matrix.i.w += x;
        this.matrix.j.w += y;
        this.matrix.k.w += z;
        this.alignTriangles();
    }
    /**
     *
     * @param {number} s
     */
    scale(s) {
        this.matrix = Matrix.scale(this.matrix, s);
        this.alignTriangles();
    }
    /**
     *
     * @param {Ray} ray
     * @returns {[number, Triangle]}
     */
    trace(ray) {
        let te = Infinity;
        let triangle = undefined;
        for (let i = 0; i < this.triangles.length; i++) {
            let t = Triangle.triangleRayIntersection(ray, this.triangles[i]);
            if (!t || t <= 0) continue;
            if (t < te) {
                te = t;
                triangle = this.triangles[i];
            }
        }
        return [te, triangle];
    }
}

const cube = new Mesh("cube.json");
//const cube2 = new Mesh("cube.json");
let film = new Film(new Vector3D(0, 0, 0), 300, 300, ctx, 2);
let camera = new Camera(film, 400);

function animate() {
    console.time("frame");
    requestAnimationFrame(animate);
    // translation every frame
    cube.translate(0, 0, 0);
    // rotation every frame
    cube.rotateXYZ(
        (-3.2 / 180) * Math.PI,
        (1 / 180) * Math.PI,
        (2 / 180) * Math.PI
    );
    // cube.scale(1.001);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (let i = 0; i < film.pixelsInSpace.length; i++) {
        let ray = camera.trace(i);
        let t = Infinity;
        meshes.forEach((m) => {
            let r = m.trace(ray);
            let te = r[0];
            let tr = r[1];
            if (te < t) {
                te = t;
                ray.color = tr.color;
            }
        });
        film.drawPixelbyIndex(i, ray.color);
    }
    console.timeEnd("frame");
}
animate();
