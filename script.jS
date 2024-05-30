//Constants
const gridHeight = 600
const gridWidth = 600
const dt = 0.1
const viscosity = 0.1
const rho = 1.0

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

let fx = createArr(gridHeight, gridWidth, 0)
let fy = createArr(gridHeight, gridWidth, 0)

for (let j = gridWidth / 2 - 50; j < gridWidth / 2 + 50; j ++) {
    for (let i = gridHeight / 2 - 50; i < gridHeight / 2 + 50; i ++) {
        fy[i][j] = 5;
        fx[i][j] = 5;
    }
}

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
            let y = j - dt * v[i][j];
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

    let h = 1 / gridHeight;

    for (let i = 1; i < gridHeight-1; i ++) {
        for (let j = 1; j < gridWidth - 1; j ++) {
            divergence[i][j] = 0.5 * h * (u[i][j+1] - u[i][j-1] + v[i+1][j] - v[i-1][j]);
            p[i][j] = 0;
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
    let h = 1 / gridHeight;
    for (let i = 1; i < gridHeight - 1; i ++) {
        for (let j = 1; j < gridWidth - 1; j ++) {
            u[i][j] -= 0.5 * h * (p[i][j+1] - p[i][j-1]);
            v[i][j] -= 0.5 * h * (p[i+1][j] - p[i-1][j]);
        }
    }
}

let canvas = document.getElementById("grid");
let ctx = canvas.getContext("2d");

function display(u, v, p) {
    let maxMagnitude = 0;
    for (let i = 0; i < gridHeight; i++) {
        for (let j = 0; j < gridWidth; j++) {
            let magnitude = Math.sqrt(u[i][j] ** 2 + v[i][j] ** 2);
            if (magnitude > maxMagnitude) {
                maxMagnitude = magnitude;
            }
        }
    }
    let normalizedMagnitudes = createArr(gridHeight, gridWidth, 0);
    for (let i = 1; i < gridHeight - 1; i ++) {
        for (let j = 1; j < gridWidth - 1; j ++) {
            let magnitude = Math.sqrt(u[i][j] ** 2 + v[i][j] ** 2);
            normalizedMagnitudes[i][j] = magnitude / maxMagnitude;
        }
    }
    let colors = getColors(normalizedMagnitudes);
    for (let i = 1; i < gridHeight; i ++) {
        for (let j = 1; j < gridWidth; j ++) {
            ctx.fillStyle = "rgb(" + colors[i][j][0] + ", 0, " + colors[i][j][1] + ")";
            ctx.fillRect(j * canvas.width / gridWidth, i * canvas.height / gridHeight, 1, 1);
        }
    }
}

function getColors(normalizedMagnitudes) {
    let colors = createArr(gridHeight, gridWidth, 0);
    for (let i = 0; i < gridHeight; i ++) {
        colors[i] = [];
        for (let j = 0; j < gridWidth; j ++) {
            colors[i][j] = [255 * normalizedMagnitudes[i][j], 255 * (1 - normalizedMagnitudes[i][j])];
        }
    }
    return colors;
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
    setInterval(simulate, 30);
}