struct Kernel {

    std::vector<float> W;

    std::vector<Vector> dWdx;
};

void compute_kernel(
        RunTime runTime,
        const float h,
        const ParticleNeighbours& particle_neighbours,
        Kernel& kernel
        )
{
    runTime.info_begin() << "Computing Kernel";
    const size_t size {particle_neighbours.origId.size()};
    const float invh = 1/h;
    const float W_fak2  = 21. / (256. * M_PI * h * h * h);
    const float dW_fak2 = 21. / (256. * M_PI * h * h * h * h);

    for (size_t pid = 0; pid < size; pid++) {

        const float len = std::sqrt(particle_neighbours.squared_length[pid]);

        // std::cout
        //     << "[DEBUG] len " << len
        //     << "pid " << pid
        //     << std::endl;

        const float q {len * invh};

        // if (q > 2.) {
        //     std::cout << "[DEBUG] outside kernel radius" << std::endl;
        //     kernel.W[pid] = 0.0;
        //     kernel.dWdx[pid] = Vector {0.0, 0.0, 0.0};
        //     continue;
        // }

        const float q3 = (q - 2.);
        const float qfac2 = q3*q3;
        const float qfac4 = qfac2*qfac2;

        float q2 = 2.*q;
        q2 += 1.;

        kernel.W[pid] = qfac4*q2*W_fak2;

        const float prefact = 10. * qfac2 * q * dW_fak2;
        kernel.dWdx[pid] = particle_neighbours.normalised_distances[pid]*prefact;
    }

    runTime.info_end();
};

void compute_nu(
        RunTime& runTime,
        const ParticleNeighbours& particle_neighbours,
        const std::vector<Vector> &velocities,
        const std::vector<float> &rho,
        const std::vector<float> &nu,
        const float h
        ) {
    // return cp.div(
    //         cp.filter(cp.scalar_product(dist, vel_ab), h, 0.0000001),
    //         cp.scalar_product(dist, dist))

    runTime.info_begin() << "Computing viscosity";
    std::vector<Vector> velocity_ab = particle_neighbours.subtract(velocities);
    auto ret = divides(
        multiply(velocity_ab, particle_neighbours.distances),
        particle_neighbours.squared_length // TODO fix this
        );

    runTime.info_end();
}

void compute_density(
        RunTime& runTime,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        std::vector<float>& densities
        ) {
    runTime.info_begin() << "Computing density";
    particle_neighbours.sum(kernel.W, densities);
    runTime.info_end();
}


void compute_du(
        RunTime& runTime,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        std::vector<Vector>& du,
        std::vector<Vector>& velocities,
        std::vector<float>& densities,
        std::vector<float>& nu,
        const float dx
        )
{
        std::vector<float> p (velocities.size(), 0);

        std::vector<float> pdd = particle_neighbours.subtract(divides(p, sqr(densities)));
        std::vector<float> temp = sum(pdd, nu);

        particle_neighbours.sum(multiply(temp, kernel.dWdx), du);
}


void compute_u(
        RunTime& runTime,
        const ParticleNeighbours& particle_neighbours,
        const Kernel & kernel,
        std::vector<Vector>& du,
        std::vector<Vector>& velocities,
        std::vector<float>& densities,
        std::vector<float>& nu,
        const float dx,
        const double dt
        )
{
  velocities = add(velocities, multiply(dt,du));
}



// def du(vel, dist, press, dens, kern, g, nb, h):
//     # monaghan eq 3.3 (P_b(rho_b)^2 + P_a(rho_a)^2
//     pdd = cp.op_ab(press/(dens**2), nb, "sub")
//     t = cp.add(pdd, nu(dist, vel, dens, nb, h))
//     return cp.sum(cp.mult_scalar_vec(t, kern), nb)

