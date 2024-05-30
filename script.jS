// Constants
const gridHeight = 300;
const gridWidth = 300;
const dt = 0.5;
const viscosity = 0.0001;
const rho = 1.0;

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

// Initializing the grid values
let u = createArr(gridHeight, gridWidth, 2);
let v = createArr(gridHeight, gridWidth, 0);
let p = createArr(gridHeight, gridWidth, 0);
let s = createArr(gridHeight, gridWidth, 0);

let fx = createArr(gridHeight, gridWidth, 0);
let fy = createArr(gridHeight, gridWidth, 0);

for (let i = gridHeight / 2 - gridHeight / 20; i < gridHeight / 2 + gridHeight / 20; i++) {
    fx[i][1] = 5;
    for (let j = 0; j < 10; j ++) {
        s[i][j] = 1.0;
    }
}

// Apply the forces to each component of the velocity vectors
function applyForces(u, v, fx, fy, dt) {
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            u[i][j] += fx[i][j] * dt;
            v[i][j] += fy[i][j] * dt;
        }
    }
}

// Advection function for a single field using the velocity vectors
function advect(u, v, field, dt) {
    let temp = createArr(gridHeight, gridWidth, 0);
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
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
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            field[i][j] = temp[i][j];
        }
    }
}

// Effect of viscosity on the movement of the fluid, reflected in the velocity field
function diffuse(u, v, viscosity, dt) {
    let tempU = createArr(gridHeight, gridWidth, 0);
    let tempV = createArr(gridHeight, gridWidth, 0);

    let diffusionCoefficient = dt * viscosity * gridHeight * gridWidth;

    // Jacobi iteration to solve the diffusion PDE
    for (let k = 0; k < 20; k++) {
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                tempU[i][j] = (u[i][j] + diffusionCoefficient * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1])) / (1 + 4 * diffusionCoefficient);
                tempV[i][j] = (v[i][j] + diffusionCoefficient * (v[i + 1][j] + v[i - 1][j] + v[i][j + 1] + v[i][j - 1])) / (1 + 4 * diffusionCoefficient);
            }
        }
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                u[i][j] = tempU[i][j];
                v[i][j] = tempV[i][j];
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
            divergence[i][j] = 0.5 * h * (u[i][j + 1] - u[i][j - 1] + v[i + 1][j] - v[i - 1][j]);
            p[i][j] = 0;
        }
    }

    // Jacobi iteration using k = 20
    for (let k = 0; k < 20; k++) {
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                temp[i][j] = 0.25 * (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - divergence[i][j]);
            }
        }
        for (let i = 1; i < gridHeight - 1; i++) {
            for (let j = 1; j < gridWidth - 1; j++) {
                p[i][j] = temp[i][j];
            }
        }
    }
}

// Project the velocity field by eliminating divergence such that the "Conservation of Mass" Navier-Stokes equation is satisfied
function projectVelocity(u, v, p) {
    let h = 1 / gridHeight;
    for (let i = 1; i < gridHeight - 1; i++) {
        for (let j = 1; j < gridWidth - 1; j++) {
            u[i][j] -= 0.5 * h * (p[i][j + 1] - p[i][j - 1]);
            v[i][j] -= 0.5 * h * (p[i + 1][j] - p[i - 1][j]);
        }
    }
}

let canvas = document.getElementById("grid");
let ctx = canvas.getContext("2d");

function display(f) {
    let field;
    if (f == "pressure") {
        field = p;
    } else if (f == "smoke") {
        field = s;
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
    let colors = getColors(normalizedMagnitudes);
    for (let i = 1; i < gridHeight; i++) {
        for (let j = 1; j < gridWidth; j++) {
            ctx.fillStyle = "rgb(" + colors[i][j][0] + ", 0, " + colors[i][j][1] + ")";
            ctx.fillRect(j * canvas.width / gridWidth, i * canvas.height / gridHeight, canvas.width / gridWidth, canvas.height / gridHeight);
        }
    }
}

function getColors(normalizedMagnitudes) {
    let colors = createArr(gridHeight, gridWidth, 0);
    for (let i = 0; i < gridHeight; i++) {
        colors[i] = [];
        for (let j = 0; j < gridWidth; j++) {
            colors[i][j] = [255 * normalizedMagnitudes[i][j], 255 * (1 - normalizedMagnitudes[i][j])];
        }
    }
    return colors;
}

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
    applyForces(u, v, fx, fy, dt);
    advect(u, v, u, dt);
    advect(u, v, v, dt);
    advect(u, v, s, dt);
    diffuse(u, v, viscosity, dt);
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