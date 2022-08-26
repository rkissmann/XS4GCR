// Copyright (c) 2017 Carmelo Evoli - MIT License

#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <sstream>

#include "XS4GCR/pid.h"
#include "XS4GCR/xs4gcr.h"

void write_table(const std::shared_ptr<XS4GCR::Spallation>& sigma,
                 const std::shared_ptr<XS4GCR::CosmicRayChart>& chart,
                 const std::string& filename_base,
                 bool doGhosts) {

    const double E_min = 10. * cgs::MeV;
    const double E_max = 1000. * cgs::TeV;
    const size_t E_size = 127;
    const XS4GCR::log_axis<double> E(E_min, E_max, E_size);
 
    std::ostringstream fname_stream;
    fname_stream << filename_base;
    if(doGhosts) {
        fname_stream << "_cumulative_";
    } else {
        fname_stream << "_direct_";
    }
    fname_stream << "Emin" << E.at(0)/ cgs::GeV;
    fname_stream << "_Emax" << E.at(E_size-1)/ cgs::GeV;
    fname_stream << "_nE" << E_size;
    fname_stream << ".txt";
    std::string filename = fname_stream.str();

    std::fstream txtfile;
    txtfile.open(filename, std::ios_base::out);
    txtfile.precision(3);

    txtfile << "0 0 0 0 ";
    for (size_t i = 0; i < E_size; ++i) {
        txtfile << std::scientific << E.at(i) / cgs::GeV;
        if (i < E_size - 1) txtfile << " ";
    }
    txtfile << "\n";

    auto particles = chart->get_particle_list();

    size_t counter = 0;
    for (auto& projectile : particles) {
        for (auto& fragment : particles) {
            if (projectile.get_A() > fragment.get_A()) {
                counter++;
                txtfile << fragment.get_Z() << " " << fragment.get_A() << " ";
                txtfile << projectile.get_Z() << " " << projectile.get_A() << " ";
                for (size_t i = 0; i < E_size; ++i) {
                    XS4GCR::channel ch(projectile, fragment);
                    txtfile << sigma->get(ch, XS4GCR::H_ISM, E.at(i), doGhosts) / cgs::mbarn;
                    if (i < E_size - 1) txtfile << " ";
                }
                txtfile << "\n";
            }
        }
    }
    std::cout << " - computed " << counter << " reactions\n";
}

int main() {
    XS4GCR::XSECS xsec;

    xsec.set_secondary_nuclei("Evoli2019");
    auto x_in = xsec.create_secondary_nuclei();

    auto crchart = xsec.create_decay_chart();

    write_table(x_in, crchart, "fragmentation_Evoli2019", false);
    write_table(x_in, crchart, "fragmentation_Evoli2019", true);
}
