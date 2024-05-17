#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>

// Сетка
const int Nx = 100;
const int Ny = 100;
const int Nz = 100;

// Источник тепла
const int hx = Nx / 2;
const int hy = Ny / 2;
const int hz = Nz / 2;
const double heat_source_value = 100.0;

const int iterations = 100;

int q = 2;

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
            std::vector<double> rbuf_send(sy * sz * q, 0);
            std::vector<double> rbuf_recv(sy * sz * q, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        rbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[sx - x][y][z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        grid[sx + q - x][y][z] = rbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                    }
                }
            }
        } else if (cx == nx - 1) {
            std::vector<double> lbuf_send(sy * sz * q, 0);
            std::vector<double> lbuf_recv(sy * sz * q, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        lbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[q - x][y][z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        grid[q - 1 - x][y][z] = lbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                    }
                }
            }
        } else if (cx % 2 == 0) {
            std::vector<double> rbuf_send(sy * sz * q, 0);
            std::vector<double> rbuf_recv(sy * sz * q, 0);
            std::vector<double> lbuf_send(sy * sz * q, 0);
            std::vector<double> lbuf_recv(sy * sz * q, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        rbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[sx - x][y][z];
                        lbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[q - x][y][z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        grid[sx + q - x][y][z] = rbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                        grid[q - 1 - x][y][z] = lbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                    }
                }
            }
        } else {
            std::vector<double> rbuf_send(sy * sz * q, 0);
            std::vector<double> rbuf_recv(sy * sz * q, 0);
            std::vector<double> lbuf_send(sy * sz * q, 0);
            std::vector<double> lbuf_recv(sy * sz * q, 0);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        rbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[sx - x][y][z];
                        lbuf_send[(y - 1) * sz * q + (z - 1) * q + x] = grid[q - x][y][z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sy * sz * q, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sy * sz * q, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
            for (int y = 1; y <= sy; y++) {
                for (int z = 1; z <= sz; z++) {
                    for (int x = 0; x < q; x++) {
                        grid[sx + q - x][y][z] = rbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                        grid[q - 1 - x][y][z] = lbuf_recv[(y - 1) * sz * q + (z - 1) * q + x];
                    }
                }
            }
        }
    }
    //y
    if (ny > 1) {
        if (cy == 0) {
            std::vector<double> rbuf_send(sx * sz * q, 0);
            std::vector<double> rbuf_recv(sx * sz * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        rbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][sy - y][z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        grid[x][sy + q - y][z] = rbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                    }
                }
            }
        } else if (cy == ny - 1) {
            std::vector<double> lbuf_send(sx * sz * q, 0);
            std::vector<double> lbuf_recv(sx * sz * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        lbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][q - y][z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        grid[x][q - 1 - y][z] = lbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                    }
                }
            }
        } else if (cx % 2 == 0) {
            std::vector<double> rbuf_send(sx * sz * q, 0);
            std::vector<double> rbuf_recv(sx * sz * q, 0);
            std::vector<double> lbuf_send(sx * sz * q, 0);
            std::vector<double> lbuf_recv(sx * sz * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        rbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][sy - y][z];
                        lbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][q - y][z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        grid[x][sy + q - y][z] = rbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                        grid[x][q - 1 - y][z] = lbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                    }
                }
            }
        } else {
            std::vector<double> rbuf_send(sx * sz * q, 0);
            std::vector<double> rbuf_recv(sx * sz * q, 0);
            std::vector<double> lbuf_send(sx * sz * q, 0);
            std::vector<double> lbuf_recv(sx * sz * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        rbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][sy - y][z];
                        lbuf_send[(x - 1) * sz * q + (z - 1) * q + y] = grid[x][q - y][z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sz * q, MPI_DOUBLE, rank - nx, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sz * q, MPI_DOUBLE, rank + nx, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int z = 1; z <= sz; z++) {
                    for (int y = 0; y < q; y++) {
                        grid[x][sy + q - y][z] = rbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                        grid[x][q - 1 - y][z] = lbuf_recv[(x - 1) * sz * q + (z - 1) * q + y];
                    }
                }
            }
        }
    }
    //z
    if (nz > 1) {
        if (cz == 0) {
            std::vector<double> rbuf_send(sx * sy * q, 0);
            std::vector<double> rbuf_recv(sx * sy * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        rbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][sz - z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        grid[x][y][sz + q - z] = rbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                    }
                }
            }
        } else if (cz == nz - 1) {
            std::vector<double> lbuf_send(sx * sy * q, 0);
            std::vector<double> lbuf_recv(sx * sy * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 1; z < q; z++) {
                        lbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][q - z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 1; z < q; z++) {
                        grid[x][y][q - 1 - z] = lbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                    }
                }
            }
        } else if (cz % 2 == 0) {
            std::vector<double> rbuf_send(sx * sy * q, 0);
            std::vector<double> rbuf_recv(sx * sy * q, 0);
            std::vector<double> lbuf_send(sx * sy * q, 0);
            std::vector<double> lbuf_recv(sx * sy * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        rbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][sz - z];
                        lbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][q - z];
                    }
                }
            }
            MPI_Send(&rbuf_send[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&lbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        grid[x][y][sz + q - z] = rbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                        grid[x][y][q - 1 - z] = lbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                    }
                }
            }
        } else {
            std::vector<double> rbuf_send(sx * sy * q, 0);
            std::vector<double> rbuf_recv(sx * sy * q, 0);
            std::vector<double> lbuf_send(sx * sy * q, 0);
            std::vector<double> lbuf_recv(sx * sy * q, 0);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        rbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][sz - z];
                        lbuf_send[(x - 1) * sy * q + (y - 1) * q + z] = grid[x][y][q - z];
                    }
                }
            }
            MPI_Recv(&lbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD, &status);
            MPI_Send(&lbuf_send[0], sx * sy * q, MPI_DOUBLE, rank - nx*ny, 1, MPI_COMM_WORLD);
            MPI_Send(&rbuf_send[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD);
            MPI_Recv(&rbuf_recv[0], sx * sy * q, MPI_DOUBLE, rank + nx*ny, 1, MPI_COMM_WORLD, &status);
            for (int x = 1; x <= sx; x++) {
                for (int y = 1; y <= sy; y++) {
                    for (int z = 0; z < q; z++) {
                        grid[x][y][sz + q - z] = rbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                        grid[x][y][q - 1 - z] = lbuf_recv[(x - 1) * sy * q + (y - 1) * q + z];
                    }
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

    int dx = (cx < Nx % nx) ? cx * sx - q: Nx % nx + cx * sx - q;
    int dy = (cy < Ny % ny) ? cy * sy - q: Ny % ny + cy * sy - q;
    int dz = (cz < Nz % nz) ? cz * sz - q: Nz % nz + cz * sz - q;

    std::vector<std::vector<std::vector<double>>> grid(sx + 2*q, std::vector<std::vector<double>>(sy + 2*q, std::vector<double>(sz + 2*q, 0.0)));
    std::vector<std::vector<std::vector<double>>> new_grid(sx + 2*q, std::vector<std::vector<double>>(sy + 2*q, std::vector<double>(sz + 2*q, 0.0)));

    initialize_grid(grid, dx, dy, dz);

    MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime();

    for (int iter = 0; iter < iterations / q; iter++) {
        exchange_boundaries(grid, rank, nx, ny, nz, sx, sy, sz, cx, cy, cz);

        for (int i = 0; i < q; i++) {
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
        }

        grid.swap(new_grid);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime() - t;

    /*
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
    */

    if (rank == 0) {
        std::cout << size << " " << t << std::endl;
    }

    MPI_Finalize();
    return 0;
}
