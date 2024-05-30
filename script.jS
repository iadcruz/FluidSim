//Constants
const gridHeight = 600
const gridWidth = 600
const dt = 0.01
const viscosity = 0.1
const rho = 1.0

const CLF = 1.9

//Function creates 2D arrays that map out the simulation grid
function createArr(x, y, val) {
    let arr = [];
    for (let i = 0; i < x; i ++) {
        arr[i] = [];
        for(let j = 0; j < y; j ++) {
            arr[i][j] = val;
        }
    }
    return arr;
}

//Initializing the grid values
let u = createArr(gridHeight, gridWidth, 0)
let v = createArr(gridHeight, gridWidth, 0)
let p = createArr(gridHeight, gridWidth, 0)

let fx = createArr(gridHeight, gridWidth, 5)
let fy = createArr(gridHeight, gridWidth, -5)

/*
for (let i = gridHeight / 2 - 25; i < gridHeight / 2 + 25; i ++) {
    fx[1][i] = 20;
} */

//Apply the forces to each component of the velocity vectors
function applyForces(u, v, fx, fy, dt) {
    for (let i = 1; i < gridHeight-1; i ++) {
        for (let j = 1; j < gridWidth-1; j ++) {
            u[i][j] += fx[i][j] * dt;
            v[i][j] += fy[i][j] * dt;
        }
    }
}

//Advects the inputted field using the velocity vectors provided
function advection(u, v, field, dt) {
    let temp = createArr(gridHeight, gridWidth, 0);
    for (let i = 1; i < gridHeight-1; i ++) {
        for (let j = 1; j < gridWidth-1; j ++) {
            let x = i - dt * u[i][j];
            let y = i - dt * v[i][j];
            x = Math.max(x, 0.5);
            x = Math.min(x, gridWidth - 1.5);
            y = Math.max(y, 0.5);
            y = Math.min(y, gridHeight - 1.5);
            let x0 = Math.floor(x);
            let y0 = Math.floor(y);
            temp[i][j] = (1-(x-x0)) * ((y-y0) * field[x0][y0+1] + (1-(y-y0)) * field[x0][y0]) + (x-x0) * ((y-y0) * field[x0+1][y0+1] + (1-(y-y0)) * field[x0+1][y0]);
        }
    }
    for (let i = 0; i < gridHeight; i ++) {
        for (let j = 0; j < gridWidth; j ++) {
            field[i][j] = temp[i][j];
        }
    }
}

//Effect of viscosity on the movement of the fluid, reflected in the velocity field
function diffuse(u, v, viscosity, dt) {
    let tempU = createArr(gridHeight, gridWidth, 0);
    let tempV = createArr(gridHeight, gridWidth, 0);

    let diffusionCoefficient = dt * viscosity * gridHeight * gridWidth;

    //Jacobi iteration to solve the diffusion PDE
    for (let k = 0; k < 5; k ++) {
        for (let i = 1; i < gridHeight - 1; i ++) {
            for (let j = 1; j < gridWidth - 1; j ++) {
                tempU[i][j] = (u[i][j] + diffusionCoefficient * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1])) / (1 + 4 * diffusionCoefficient);
                tempV[i][j] = (v[i][j] + diffusionCoefficient * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1])) / (1 + 4 * diffusionCoefficient);
            }
        }
        for (let i = 1; i < gridHeight - 1; i ++) {
            for (let j = 1; j < gridWidth - 1; j ++) {
                u[i][j] = tempU[i][j];
                v[i][j] = tempV[i][j];
            }
        }
    }
}

//Pressure calculation using the Jacobi iterative method to solve the Poisson equation
function pressureCalculation(p, u, v, rho, dt) {
    let divergence = createArr(gridHeight, gridWidth, 0);
    let temp = createArr(gridHeight, gridWidth, 0);

    let h = CLF * 1 / gridHeight;

    for (let i = 1; i < gridHeight-1; i ++) {
        for (let j = 1; j < gridWidth - 1; j ++) {
            divergence[i][j] = 0.5 * h * (u[i+1][j] - u[i-1][j] + v[i][j+1] - v[i][j-1])
            p[i][j] = 0
        }
    }

    //Jacobi iteration using k = 5
    for (let k = 0; k < 5; k ++) {
        for (let i = 1; i < gridHeight - 1; i ++) {
            for (let j = 1; j < gridWidth - 1; j ++) {
                temp[i][j] = 0.25 * (p[i+1][j] + p[i-1][j] + p[i][j+1] + p[i][j-1] - divergence[i][j]);
            }
        }
        for (let i = 1; i < gridHeight - 1; i ++) {
            for (let j = 1; j < gridWidth - 1; j ++) {
                p[i][j] = temp[i][j];
            }
        }
    }
}

//Project the velocity field by eliminating divergence such that the "Conservation of Mass" Navier-Stokes equation is satisfied
function projectVelocity(u, v, p) {
    let h = CLF * 1 / gridHeight;
    for (let i = 1; i < gridHeight - 1; i ++) {
        for (let j = 1; j < gridWidth - 1; j ++) {
            u[i][j] -= 0.5 * h * (p[i+1][j] - p[i-1][j]);
            v[i][j] -= 0.5 * h * (p[i][j+1] - p[i][j-1]);
        }
    }
}

let canvas = document.getElementById("grid");
let ctx = canvas.getContext("2d");

function display(u, v, p) {

    const colors = ['#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000'];

    let maxMagnitude = 0;
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            const magnitude = Math.sqrt(u[i][j] ** 2 + v[i][j] ** 2);
            maxMagnitude = Math.max(maxMagnitude, magnitude);
        }
    }
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw velocity vectors with color representing magnitude
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            const magnitude = Math.sqrt(u[i][j] ** 2 + v[i][j] ** 2);
            const normalizedMagnitude = magnitude / maxMagnitude;
            const colorIndex = Math.floor(normalizedMagnitude * (colors.length - 1));
            const colorFraction = (normalizedMagnitude * (colors.length - 1)) % 1;
            const color1 = colors[colorIndex];
            const color2 = colors[colorIndex + 1] || colors[colorIndex];
            const color = interpolateColor(color1, color2, colorFraction);
            ctx.fillStyle = color;
            ctx.beginPath();
            ctx.arc(j, i, 2, 0, 2 * Math.PI);
            ctx.fill();
        }
    }
}

function interpolateColor(color1, color2, fraction) {
    const r1 = parseInt(color1.substring(1, 3), 16);
    const g1 = parseInt(color1.substring(3, 5), 16);
    const b1 = parseInt(color1.substring(5, 7), 16);
    const r2 = parseInt(color2.substring(1, 3), 16);
    const g2 = parseInt(color2.substring(3, 5), 16);
    const b2 = parseInt(color2.substring(5, 7), 16);

    const r = Math.round(r1 + (r2 - r1) * fraction);
    const g = Math.round(g1 + (g2 - g1) * fraction);
    const b = Math.round(b1 + (b2 - b1) * fraction);

    return '#' + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

function componentToHex(c) {
    const hex = c.toString(16);
    return hex.length == 1 ? '0' + hex : hex;
}

//Calls the individual functions to run the simulation
function simulate() {
    applyForces(u, v, fx, fy, dt);
    advection(u, v, u, dt);
    advection(u, v, v, dt);
    diffuse(u, v, viscosity, dt);
    pressureCalculation(p, u, v, rho, dt);
    projectVelocity(u, v, p);
    display(u, v, p);
}

function startSimulation(n) {
    setInterval(simulate, 100)
}