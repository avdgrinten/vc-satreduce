#include <bnb.hpp>
#include <chrono>
#include <fstream>
#include <graph.hpp>
#include <input-reader.hpp>
#include <iostream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    po::options_description desc("options");
    // clang-format off
    desc.add_options()
        // Options that control output.
        ("help,h", "Display help message")
        ("verbose,v", "Debug output")
        ("stats", "Display statistics")
        ("print-solution", "Print the solution")

        // Options that control solver behavior.
        ("no-initial-solution", "Do not compute heuristic initial solution")
        ("lp-reduction", "Use LP reduction")
        ("packing", po::value<bool>()->default_value(true, "true"),
            "Use packing")
        ("funnel", po::value<bool>()->default_value(false, "false"),
            "Use funnel in all branches of the search tree")
        ("sat-conflict-limit,", po::value<int>(), "Conflict limit for SAT solver")
        ("sat-invalidation-count", po::value<int>(), "")
        ("sat-factor", po::value<double>(), "")
        ("no-sat,", "Disable SAT solver")
        ("vertex-encoding,", "")
        ("constraint-totalizer,", "")
        ("only-sat", "")
        ("with-reductions", "")
        ("ub-to-lb", "")
        ("check-fixed-yc", "")
        ;
    // clang-format on

    po::options_description full_desc("all options");
    // clang-format off
    full_desc.add_options()
        ("instance", po::value<std::string>(), "Specify the instance file")
        ;
    // clang-format on
    full_desc.add(desc);

    po::positional_options_description pdesc;
    pdesc.add("instance", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(full_desc).positional(pdesc).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "usage: vc-bnb [options...] instance" << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    std::ifstream stream{vm["instance"].as<std::string>()};
    Input_Reader reader;
    Graph graph = reader.read_instance(stream);
    vc_bnb::Solver sol(graph);

    // init flags
    sol.flags = Flags(vm);

    auto start = std::chrono::high_resolution_clock::now();

    bool success = sol.solve();

    if (vm.count("stats")) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        std::cout << "total:\n";
        std::cout << "    graph_size: " << sol.get_start_vertex_count() << "\n";
        std::cout << "    total_time: " << diff.count() << "\n";
        std::cout << "    total_branches: " << sol.iteration_count << "\n";
        if (success)
            std::cout << "    vertex_cover: " << sol.best_solution.size() << std::endl;
        sol.print_stats();
    }

    if (success && vm.count("print-solution")) {
        std::cout << "vc:\n";
        int n = 1;
        for (auto v : sol.best_solution.picked_vertices) {
            if (v)
                std::cout << "  - " << n << "\n";

            n++;
        }
        std::cout << std::flush;
    }
}
