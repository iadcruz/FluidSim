// Constants
const gridHeight = 200;
const gridWidth = 200;
const dt = 0.5;
const viscosity = 0.1;
const rho = 0.1;

// Function creates 2D arrays that map out the simulation grid
function createArr(x, y, val) {
    let arr = [];
    for (let i = 0; i < x; i++) {
        arr[i] = [];
        for (let j = 0; j < y; j++) {
            arr[i][j] = val;
        }
    }
    return arr;
}

//Creating Canvas
let canvas = document.getElementById("grid");
let ctx = canvas.getContext("2d");

// Initializing the grid values
let u = createArr(gridHeight, gridWidth, 0);
let v = createArr(gridHeight, gridWidth, 0);
let p = createArr(gridHeight, gridWidth, 0);
let s = createArr(gridHeight, gridWidth, 0);
let d = createArr(gridHeight, gridWidth, 0);
let fx = createArr(gridHeight, gridWidth, 0);
let fy = createArr(gridHeight, gridWidth, 1);
let obstacleGrid = createArr(gridHeight, gridWidth, false);

// Set initial conditions
for (let i = gridHeight / 2 - 30; i < gridHeight / 2 + 30; i++) {
    for (let j = 1; j < gridWidth / 10; j++) {
        fy[i][j] = 10;
        u[i][j] = 20;
        v[i][j] = 40;
        s[i][j] = 1.0;
    }
}

//Create obstacle

//Rectangle

function rectangle() {
    for (let i = 75; i < 125; i ++) {
        for (let j = 75; j < 125; j ++) {
            obstacleGrid[i][j] = true;
        }
    }
}

let time = 0;
const MaxDistanceDisplaced = 30;

function moveobstacle(obstacle) {
    time += 1;
    let tempamplitude = Math.floor(MaxDistanceDisplaced * Math.sin(time));
    obstacleGrid = createArr(gridHeight, gridWidth, false);
    if (obstacle == "rectangle") {
        for (let i = 75 + tempamplitude; i < 125 + tempamplitude; i ++) {
            for (let j = 75; j < 125; j ++) {
                obstacleGrid[i][j] = true;
            }
        }
    }
}

//Circle
let r = 25;
let yC = 100;
let xC = 100;

/* for (let j = yC - r; j < yC + r; j ++) {
    xL = Math.floor(Math.sqrt(Math.pow(r, 2) - Math.pow(j - yC, 2)));
    for (let i = xC - xL; i < xC + xL; i ++) {
        obstacleGrid[i][j] = true;
    }
}
*/

// Apply the forces to each component of the velocity vectors
function applyForces(u, v, fx, fy, dt) {
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            if (!obstacleGrid[i][j]) {
                u[i][j] += fx[i][j] * dt;
                v[i][j] += fy[i][j] * dt;
            }
        }
    }
}

// Updates smoke densities based on sources
function sourceSmoke(u, v, d, s, dt) {
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            if (!obstacleGrid[i][j]) {
                d[i][j] += s[i][j] * dt;
            }
        }
    }
}

// Advection function for a single field using the velocity vectors
function advect(u, v, field, dt) {
    let temp = createArr(gridHeight, gridWidth, 0);
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            if (!obstacleGrid[i][j]) {
                let x = i - dt * u[i][j];
                let y = j - dt * v[i][j];
                x = Math.max(x, 0.5);
                x = Math.min(x, gridWidth - 1.5);
                y = Math.max(y, 0.5);
                y = Math.min(y, gridHeight - 1.5);
                let x0 = Math.floor(x);
                let y0 = Math.floor(y);
                let x1 = x0 + 1;
                let y1 = y0 + 1;
                let s1 = x - x0;
                let t1 = y - y0;
                let s0 = 1 - s1;
                let t0 = 1 - t1;
                temp[i][j] = s0 * (t0 * field[x0][y0] + t1 * field[x0][y1]) + s1 * (t0 * field[x1][y0] + t1 * field[x1][y1]);
            }
        }
    }
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            if (!obstacleGrid[i][j]) {
                field[i][j] = temp[i][j];
            }
        }
    }
}

// Effect of viscosity on the movement of the fluid, reflected in the velocity field and smoke density
function diffuse(field, viscosity, dt) {
    let temp = createArr(gridHeight, gridWidth, 0);
    let diffusionCoefficient = dt * viscosity * gridHeight * gridWidth;

    // Jacobi iteration to solve the diffusion PDE
    for (let k = 0; k < 20; k++) {
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                if (!obstacleGrid[i][j]) {
                    temp[i][j] = (field[i][j] + diffusionCoefficient * (field[i + 1][j] + field[i - 1][j] + field[i][j + 1] + field[i][j - 1])) / (1 + 4 * diffusionCoefficient);
                }
            }
        }
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                if (!obstacleGrid[i][j]) {
                    field[i][j] = temp[i][j];
                }
            }
        }
    }
}

// Pressure calculation using the Jacobi iterative method to solve the Poisson equation
function pressureCalculation(p, u, v, rho, dt) {
    let divergence = createArr(gridHeight, gridWidth, 0);
    let temp = createArr(gridHeight, gridWidth, 0);
    let h = 1 / gridHeight;

    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            if (!obstacleGrid[i][j]) {
                divergence[i][j] = 0.5 * h * (u[i][j + 1] - u[i][j - 1] + v[i + 1][j] - v[i - 1][j]);
                p[i][j] = 0;
            }
        }
    }

    // Jacobi iteration using k = 20
    for (let k = 0; k < 20; k++) {
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                if (!obstacleGrid[i][j]) {
                    temp[i][j] = 0.25 * (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - divergence[i][j]);
                }
            }
        }
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                if (!obstacleGrid[i][j]) {
                    p[i][j] = temp[i][j];
                }
            }
        }
    }
}

// Project the velocity field by eliminating divergence such that the "Conservation of Mass" Navier-Stokes equation is satisfied
function projectVelocity(u, v, p) {
    let h = 1 / gridHeight;
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            if (!obstacleGrid[i][j]) {
                u[i][j] -= 0.5 * h * (p[i][j + 1] - p[i][j - 1]);
                v[i][j] -= 0.5 * h * (p[i + 1][j] - p[i - 1][j]);
            }
        }
    }
}

function display(f) {
    let field;
    if (f == "pressure") {
        field = p;
    } else if (f == "smoke") {
        field = d;
    } else if (f == "velocity") {
        field = calculateVelocity();
    }
    let maxMagnitude = 0;
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            let magnitude = field[i][j];
            if (magnitude > maxMagnitude) {
                maxMagnitude = magnitude;
            }
        }
    }
    let normalizedMagnitudes = createArr(gridHeight, gridWidth, 0);
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            normalizedMagnitudes[i][j] = field[i][j] / maxMagnitude;
        }
    }
    // let colors = getbasic(field);
    let colors = createArr(gridHeight, gridWidth, 0);
    for (let i = 0; i < gridHeight; i ++) {
        for (let j = 0; j < gridWidth; j ++) {
            colors[i][j] = getColor(normalizedMagnitudes[i][j]);
        }
    }
    for (let i = 1; i < gridHeight; i++) {
        for (let j = 1; j < gridWidth; j++) {
            //ctx.fillStyle = colors[i][j];
            ctx.fillStyle = "rgb(" + colors[i][j][0] + ", " + colors[i][j][1] + ", " + colors[i][j][2] + ")";
            ctx.fillRect(j * canvas.width / gridWidth, i * canvas.height / gridHeight, canvas.width / gridWidth, canvas.height / gridHeight);
        }
    }
}

function getbasic(field) {
    let colors = createArr(gridHeight, gridWidth, 0);
    for (let i = 0; i < gridHeight; i ++) {
        for (let j = 0; j < gridWidth; j ++) {
            if (field[i][j] != 0) {
                colors[i][j] = "gray";
            }
            else {
                colors[i][j] = "white";
            }
        }
    }
    return colors;
}

function interpolateColor(value, color1, color2) {
    let r = Math.floor(color1.r + value * (color2.r - color1.r));
    let g = Math.floor(color1.g + value * (color2.g - color1.g));
    let b = Math.floor(color1.b + value * (color2.b - color1.b));
    return [r, g, b];
}

function getColor(normalized) {
    let jetColors = [
        { pos: 0.0, color: { r: 0, g: 0, b: 255 } },
        { pos: 0.35, color: { r: 0, g: 255, b: 255 } },
        { pos: 0.5, color: { r: 0, g: 255, b: 0 } },
        { pos: 0.65, color: { r: 255, g: 255, b: 0 } },
        { pos: 1.0, color: { r: 255, g: 0, b: 0 } }
    ];

    for (let i = 0; i < jetColors.length - 1; i++) {
        if (normalized >= jetColors[i].pos && normalized <= jetColors[i + 1].pos) {
            let t = (normalized - jetColors[i].pos) / (jetColors[i + 1].pos - jetColors[i].pos);
            return interpolateColor(t, jetColors[i].color, jetColors[i + 1].color);
        }
    }
    return jetColors[jetColors.length - 1].color;
}

//Calculates magnitude of the velocity field when displaying "velocity"
function calculateVelocity() {
    let arr = createArr(gridHeight, gridWidth, 0);
    for (let i = 0; i < gridHeight; i++) {
        arr[i] = [];
        for (let j = 0; j < gridWidth; j++) {
            arr[i][j] = Math.sqrt(u[i][j] ** 2 + v[i][j] ** 2);
        }
    }
    return arr;
}

// Calls the individual functions to run the simulation
function simulate(field) {
    moveobstacle("rectangle");
    applyForces(u, v, fx, fy, dt);
    sourceSmoke(u, v, d, s, dt);
    advect(u, v, u, dt);
    advect(u, v, v, dt);
    advect(u, v, d, dt);
    diffuse(u, viscosity, dt);
    diffuse(v, viscosity, dt);
    diffuse(d, viscosity, dt);
    pressureCalculation(p, u, v, rho, dt);
    projectVelocity(u, v, p);
    display(field);
}

function startSimulation() {
    let field;
    if (document.getElementById("velocity").checked) {
        field = "velocity";
    } else if (document.getElementById("pressure").checked) {
        field = "pressure";
    } else if (document.getElementById("smoke").checked) {
        field = "smoke";
    }
    setInterval(() => simulate(field), 30);
}