#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

// Сетка
const int Nx = 200;
const int Ny = 200;
const int Nz = 200;

// Источник тепла
const int hx = Nx / 2;
const int hy = Ny / 2;
const int hz = Nz / 2;
const double heat_source_value = 100.0;

const int iterations = 100;

const double k = .01;

void initialize_grid(std::vector<std::vector<std::vector<double>>>& grid, int dx, int dy, int dz) {
    for (int x = 0; x < grid.size(); x++) {
        for (int y = 0; y < grid[0].size(); y++) {
            for (int z = 0; z < grid[0][0].size(); z++) {
                if (x + dx == hx && y + dy== hy && z + dz== hz) {
                    grid[x][y][z] = heat_source_value;
                } else {
                    grid[x][y][z] = 0.0;
                }
            }
        }
    }
}

// Печатает слой сетки
void print_slice(const std::vector<std::vector<std::vector<double>>>& grid, int z) {
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            std::cout << grid[x][y][z] << " ";
        }
        std::cout << std::endl;
    }
}

void exchange_boundaries(std::vector<std::vector<std::vector<double>>>& grid, int rank, int nx, int ny, int nz,
                         int sx, int sy, int sz, int cx, int cy, int cz) {
    MPI_Status status;
    //x
    if (nx > 1) {
        if (cx == 0) {
            std::vector<double> rbuf_send(sy * sz, 0);
            std::vector<double> rbuf_recv(sy * sz, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sy + y - 1] = grid[sx][y][z];
                }
            }
            MPI_Send(&rbuf_send[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    grid[sx + 1][y][z] = rbuf_recv[(z - 1) * sy + y - 1];
                }
            }
        } else if (cx == nx - 1) {
            std::vector<double> lbuf_send(sy * sz, 0);
            std::vector<double> lbuf_recv(sy * sz, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    lbuf_send[(z - 1) * sy + y - 1] = grid[1][y][z];
                }
            }
            MPI_Recv(&lbuf_recv[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    grid[0][y][z] = lbuf_recv[(z - 1) * sy + y - 1];
                }
            }
        } else if (cx % 2 == 0) {
            std::vector<double> rbuf_send(sy * sz, 0);
            std::vector<double> rbuf_recv(sy * sz, 0);
            std::vector<double> lbuf_send(sy * sz, 0);
            std::vector<double> lbuf_recv(sy * sz, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sy + y - 1] = grid[sx][y][z];
                    lbuf_send[(z - 1) * sy + y - 1] = grid[1][y][z];
                }
            }
            MPI_Send(&rbuf_send[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    grid[sx + 1][y][z] = rbuf_recv[(z - 1) * sy + y - 1];
                    grid[0][y][z] = lbuf_recv[(z - 1) * sy + y - 1];
                }
            }
        } else {
            std::vector<double> rbuf_send(sy * sz, 0);
            std::vector<double> rbuf_recv(sy * sz, 0);
            std::vector<double> lbuf_send(sy * sz, 0);
            std::vector<double> lbuf_recv(sy * sz, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sy + y - 1] = grid[sx][y][z];
                    lbuf_send[(z - 1) * sy + y - 1] = grid[1][y][z];
                }
            }
            MPI_Recv(&lbuf_recv[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    grid[sx + 1][y][z] = rbuf_recv[(z - 1) * sy + y - 1];
                    grid[0][y][z] = lbuf_recv[(z - 1) * sy + y - 1];
                }
            }
        }
    }
    //y
    if (ny > 1) {
        if (cy == 0) {
            std::vector<double> rbuf_send(sx * sz, 0);
            std::vector<double> rbuf_recv(sx * sz, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sx + x - 1] = grid[x][sy][z];
                }
            }
            MPI_Send(&rbuf_send[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    grid[x][sy + 1][z] = rbuf_recv[(z - 1) * sx + x - 1];
                }
            }
        } else if (cy == ny - 1) {
            std::vector<double> lbuf_send(sx * sz, 0);
            std::vector<double> lbuf_recv(sx * sz, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    lbuf_send[(z - 1) * sx + x - 1] = grid[x][1][z];
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    grid[x][0][z] = lbuf_recv[(z - 1) * sx + x - 1];
                }
            }
        } else if (cx % 2 == 0) {
            std::vector<double> rbuf_send(sx * sz, 0);
            std::vector<double> rbuf_recv(sx * sz, 0);
            std::vector<double> lbuf_send(sx * sz, 0);
            std::vector<double> lbuf_recv(sx * sz, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sx + x - 1] = grid[x][sy][z];
                    lbuf_send[(z - 1) * sx + x - 1] = grid[x][1][z];
                }
            }
            MPI_Send(&rbuf_send[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    grid[x][sy + 1][z] = rbuf_recv[(z - 1) * sx + x - 1];
                    grid[x][0][z] = lbuf_recv[(z - 1) * sx + x - 1];
                }
            }
        } else {
            std::vector<double> rbuf_send(sx * sz, 0);
            std::vector<double> rbuf_recv(sx * sz, 0);
            std::vector<double> lbuf_send(sx * sz, 0);
            std::vector<double> lbuf_recv(sx * sz, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    rbuf_send[(z - 1) * sx + x - 1] = grid[x][sy][z];
                    lbuf_send[(z - 1) * sx + x - 1] = grid[x][1][z];
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    grid[x][sy + 1][z] = rbuf_recv[(z - 1) * sx + x - 1];
                    grid[x][0][z] = lbuf_recv[(z - 1) * sx + x - 1];
                }
            }
        }
    }
    //z
    if (nz > 1) {
        if (cz == 0) {
            std::vector<double> rbuf_send(sx * sy, 0);
            std::vector<double> rbuf_recv(sx * sy, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    rbuf_send[(y - 1) * sx + x - 1] = grid[x][y][sz];
                }
            }
            MPI_Send(&rbuf_send[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    grid[x][y][sz + 1] = rbuf_recv[(y - 1) * sx + x - 1];
                }
            }
        } else if (cz == nz - 1) {
            std::vector<double> lbuf_send(sx * sy, 0);
            std::vector<double> lbuf_recv(sx * sy, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    lbuf_send[(y - 1) * sx + x - 1] = grid[x][y][1];
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    grid[x][y][0] = lbuf_recv[(y - 1) * sx + x - 1];
                }
            }
        } else if (cz % 2 == 0) {
            std::vector<double> rbuf_send(sx * sy, 0);
            std::vector<double> rbuf_recv(sx * sy, 0);
            std::vector<double> lbuf_send(sx * sy, 0);
            std::vector<double> lbuf_recv(sx * sy, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    rbuf_send[(y - 1) * sx + x - 1] = grid[x][y][sz];
                    lbuf_send[(y - 1) * sx + x - 1] = grid[x][y][1];
                }
            }
            MPI_Send(&rbuf_send[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    grid[x][y][sz + 1] = rbuf_recv[(y - 1) * sx + x - 1];
                    grid[x][y][0] = lbuf_recv[(y - 1) * sx + x - 1];
                }
            }
        } else {
            std::vector<double> rbuf_send(sx * sy, 0);
            std::vector<double> rbuf_recv(sx * sy, 0);
            std::vector<double> lbuf_send(sx * sy, 0);
            std::vector<double> lbuf_recv(sx * sy, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    rbuf_send[(y - 1) * sx + x - 1] = grid[x][y][sz];
                    lbuf_send[(y - 1) * sx + x - 1] = grid[x][y][1];
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    grid[x][y][sz + 1] = rbuf_recv[(y - 1) * sx + x - 1];
                    grid[x][y][0] = lbuf_recv[(y - 1) * sx + x - 1];
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc < 4) return 0;

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int nz = atoi(argv[3]);
    
    if (nx * ny * nz != size) return 0;

    int cx = rank % nx;
    int cy = (rank / nx) % ny;
    int cz = rank / (nx * ny);

    int sx = (cx < Nx % nx) ? Nx / nx + 1 : Nx / nx;
    int sy = (cy < Ny % ny) ? Ny / ny + 1 : Ny / ny;
    int sz = (cz < Nz % nz) ? Nz / nz + 1 : Nz / nz;

    int dx = (cx < Nx % nx) ? cx * sx - 1: Nx % nx + cx * sx - 1;
    int dy = (cy < Ny % ny) ? cy * sy - 1: Ny % ny + cy * sy - 1;
    int dz = (cz < Nz % nz) ? cz * sz - 1: Nz % nz + cz * sz - 1;

    std::vector<std::vector<std::vector<double>>> grid(sx + 2, std::vector<std::vector<double>>(sy + 2, std::vector<double>(sz + 2, 0.0)));
    std::vector<std::vector<std::vector<double>>> new_grid(sx + 2, std::vector<std::vector<double>>(sy + 2, std::vector<double>(sz + 2, 0.0)));

    initialize_grid(grid, dx, dy, dz);

    MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime();

    for (int iter = 0; iter < iterations; iter++) {
        exchange_boundaries(grid, rank, nx, ny, nz, sx, sy, sz, cx, cy, cz);

        for (int x = 1; x <= sx; x++) {
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    if (x + dx == hx && y + dy== hy && z + dz== hz) {
                        new_grid[x][y][z] = heat_source_value;
                        continue;
                    }
                    new_grid[x][y][z] = grid[x][y][z] + k*(grid[x+1][y][z] + grid[x-1][y][z] + grid[x][y+1][z] + 
                                                           grid[x][y-1][z] + grid[x][y][z+1] + grid[x][y][z-1] - 6 * grid[x][y][z]);
                }
            }
        }

        grid.swap(new_grid);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;

    if (rank == 0) {
        std::vector<std::vector<std::vector<double>>> full_grid(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz, 0.0)));
        std::vector<double> buf(Nx * Ny * Nz);
        for (int x = 1; x <= sx; x++) {
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    full_grid[x + dx][y + dy][z + dz] = grid[x][y][z];
                }
            }
        }
        MPI_Status status;
        for (int r = 1; r < size; r++) {
            int rcx = r % nx;
            int rcy = (r / nx) % ny;
            int rcz = r / (nx * ny);

            int rsx = (rcx < Nx % nx) ? Nx / nx + 1 : Nx / nx;
            int rsy = (rcy < Ny % ny) ? Ny / ny + 1 : Ny / ny;
            int rsz = (rcz < Nz % nz) ? Nz / nz + 1 : Nz / nz;

            int rdx = (rcx < Nx % nx) ? rcx * rsx - 1: Nx % nx + rcx * rsx - 1;
            int rdy = (rcy < Ny % ny) ? rcy * rsy - 1: Ny % ny + rcy * rsy - 1;
            int rdz = (rcz < Nz % nz) ? rcz * rsz - 1: Nz % nz + rcz * rsz - 1;
            MPI_Recv(&buf[0], rsx * rsy * rsz, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= rsx; x++) {
                for (int y = 1; y <= rsy; y++) {
                    for (int z = 1; z <= rsz; z++) {
                        full_grid[x + rdx][y + rdy][z + rdz] = buf[(x - 1) * rsy * rsz + (y - 1) * rsz + z - 1];
                    }
                }
            }
        }
        
        //print_slice(full_grid, hz);
    } else {
        std::vector<double> buf(sx * sy * sz);
        for (int x = 1; x <= sx; x++) {
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    buf[(x - 1) * sy * sz + (y - 1) * sz + z - 1] = grid[x][y][z];
                }
            }
        }
        MPI_Send(&buf[0], sx * sy * sz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        std::cout << size << " " << t << std::endl;
    }

    MPI_Finalize();
    return 0;
}
