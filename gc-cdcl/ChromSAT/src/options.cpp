#include "options.hpp"

#include <tclap/CmdLine.h>

#include <numeric>

namespace gc
{

struct argbase {
    virtual ~argbase() {}
    virtual void assign() = 0;
};

template <typename Opt, typename ClapArg, typename E = void>
struct arg : public argbase {
    ClapArg carg;
    Opt& opt;

    template <typename... T>
    arg(TCLAP::CmdLine& cmd, Opt& opt, T&&... args)
        : carg(std::forward<T>(args)...)
        , opt(opt)
    {
        cmd.add(carg);
    }

    virtual void assign() override { opt = carg.getValue(); }
};

template <typename Opt, typename ClapArg>
struct arg<Opt, ClapArg, typename std::enable_if<std::is_enum<Opt>{}>::type>
    : public argbase {
    ClapArg carg;
    Opt& opt;

    template <typename... T>
    arg(TCLAP::CmdLine& cmd, Opt& opt, T&&... args)
        : carg(std::forward<T>(args)...)
        , opt(opt)
    {
        cmd.add(carg);
    }

    virtual void assign() override
    {
        opt = static_cast<typename std::remove_reference<Opt>::type>(
            carg.getValue());
    }
};

struct cmdline {
    TCLAP::CmdLine cmd;
    std::vector<std::unique_ptr<argbase>> args;

    cmdline(const std::string& message, const char delimiter = ' ',
        const std::string& version = "none", bool helpAndVersion = true)
        : cmd(message, delimiter, version, helpAndVersion)
    {
    }

    template <typename ClapArg, typename Opt, typename... T>
    void add(Opt& opt, T&&... clapargs)
    {
        args.emplace_back(std::move(std::make_unique<arg<Opt, ClapArg>>(
            cmd, opt, std::forward<T>(clapargs)...)));
    }

    void parse(int argc, char* argv[])
    {
        cmd.parse(argc, argv);
        for (auto& arg : args)
            arg->assign();
    }
};

options parse(int argc, char* argv[])
{
    using namespace TCLAP;
    using namespace std::string_literals;
    cmdline cmd("GraphColoring", ' ');

    options opt;
    opt.cmdline = std::accumulate(argv, argv + argc, ""s,
        [&](std::string acc, const char* arg) { return acc + " " + arg; });

    cmd.add<UnlabeledValueArg<std::string>>(opt.instance_file, "file",
        "instance file", true, "data/DIMACS_cliques/brock200_1.clq", "string");
    cmd.add<UnlabeledValueArg<std::string>>(
        opt.solution_file, "solfile", "solution file", false, "", "string");
    cmd.add<SwitchArg>(opt.trace, "", "trace", "enable minicsp tracing", false);
    cmd.add<ValueArg<int>>(opt.learning, "", "learning",
        "CDCLeaning & explanation level [0-1]", false, 1, "int");
    cmd.add<SwitchArg>(
        opt.xvars, "", "xvars", "add x (color) variables to the model", false);
    cmd.add<ValueArg<int>>(
        opt.polarity, "", "polarity", "polarity policy", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.ordering, "", "ordering",
        "clique finding heuristic [0-4]", false, 3, "int");
    cmd.add<ValueArg<int>>(opt.ordering_low_degree, "", "ord-low-degree",
        "Use low degree information to improve clique ordering", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.boundalg, "", "bound",
        "lower bound algorithm [0-3]", false, 2, "int");
    cmd.add<ValueArg<int>>(opt.cliquealg, "", "cliquealg",
        "Clique heuristic to use during search", false, 0, "0-1");
    cmd.add<ValueArg<int>>(opt.cliquemargin, "", "cliquemargin",
        "Keep cliques of size lb-x", false, 0, "int>=0");
    cmd.add<SwitchArg>(opt.prune, "", "prune", "enable positive pruning", false);
    cmd.add<SwitchArg>(
        opt.enurp, "", "enurp", "enable negative pruning", false);
    cmd.add<SwitchArg>(opt.adaptive, "", "adaptive",
        "Switch between CLIQUES and declared bound policy dynamically", true);
    cmd.add<ValueArg<int>>(opt.branching, "", "branching",
        "Variable branching heuristic [0-14]", false, 1, "int");
    cmd.add<SwitchArg>(opt.brelaz_first, "", "brelaz-first",
        "Use brelaz for 1e5 conflicts before switching to chosen heuristic",
        false);
    cmd.add<SwitchArg>(opt.phase_saving, "", "phase-saving",
        "Phase saving in branching",
        false);
    cmd.add<SwitchArg>(opt.branching_low_degree, "", "branch-low-degree",
        "Use low degree information to improve branching", false);
    cmd.add<ValueArg<int>>(opt.cliquelimit, "", "cliquelimit",
        "Maximum number of cliques in the lower bound algorithm during "
        "preprocessing",
        false, 1000, "int");
    cmd.add<ValueArg<int>>(opt.strategy, "", "strategy",
        "Solution strategy [0=BNB-1=bottom-up-2=top-down-3=preprocessing "
        "only-4=top-down as in (Verma et al.)-5=color-6=test]",
        false, 0, "int");
    cmd.add<ValueArg<int>>(opt.preprocessing, "", "preprocessing",
        "Level of preprocessing [0=none-1=low-degree-2=low-degree (sparse)]", false, 1, "int");
    cmd.add<SwitchArg>(opt.dominance, "", "dominance",
        "enable neighborhood-based dominance", false);
    cmd.add<SwitchArg>(opt.indset_constraints, "", "indset",
        "reduce by independent set constraints", false);
    cmd.add<ValueArg<int>>(opt.fillin, "", "fillin",
        "compute fill-in of graph (0:no,1:prop,2:decomp", false, 0, "int");
    cmd.add<SwitchArg>(opt.indset_lb, "", "indset-lb",
        "compute IS-based lower bound", false);
    cmd.add<ValueArg<std::string>>(opt.format, "", "format",
        "File format", false, "dimacs", "string");
    cmd.add<SwitchArg>(
        opt.checksolution, "", "checksolution", "checks the coloring", false);
    cmd.add<SwitchArg>(
        opt.printsolution, "", "printsolution", "prints the coloring", false);
    cmd.add<ValueArg<int>>(opt.sdsaturiter, "", "sdsaturiter",
        "# of sparse dsatur iterations", false, 1, "int");
    cmd.add<ValueArg<int>>(opt.ddsaturiter, "", "ddsaturiter",
        "# of dense dsatur iterations", false, 1, "int");
    cmd.add<ValueArg<std::string>>(opt.convert, "", "convert",
        "output in <file> using dimacs format", false, "", "string");
    cmd.add<SwitchArg>(opt.equalities, "", "equalities",
        "add implied number of equalitiy constraints", false);
    cmd.add<SwitchArg>(
        opt.dsatur, "", "dsatur", "use dsatur instead of minicsp", false);
    cmd.add<SwitchArg>(opt.ubfocus, "", "ubfocus",
        "do not try hard to get lb's [default=false]", false);
    cmd.add<SwitchArg>(opt.maxclique, "", "maxclique",
        "use dOmega maximum clique in preprocessing [default=false]", false);
    cmd.add<ValueArg<double>>(opt.domegatime, "", "domegatime",
        "timeout for maxclique algorithm", false, -1.0, "double");
    cmd.add<ValueArg<int>>(opt.memlimit, "", "memlimit",
        "does not compute dense graph above this limit and downgrade reasoning "
        "on large instances (default: no limit)",
        false, -1, "int");
    cmd.add<ValueArg<int>>(opt.myciellimit, "", "myciellimit",
        "does not use myciel on graphs larger than the limit (default: no limit)",
        false, -1, "int");
    cmd.add<SwitchArg>(opt.triangle_up, "", "triangleup",
        "do partition-aware unit propagation", false);
    cmd.add<ValueArg<int>>(opt.samplebase, "", "samplebase",
        "size of the sampling base for max clique (default = 512)", false, 512,
        "int");
    cmd.add<ValueArg<int>>(opt.probewidth, "", "probewidth",
        "size of the max probe width (default = 64)", false, 64, "int");
    cmd.add<ValueArg<int>>(
        opt.core, "", "core", "Core type []", false, 3, "int");
    cmd.add<ValueArg<int>>(opt.idsaturlimit, "", "idsaturlimit",
        "solve limit in idsatur", false, -1, "int");
    cmd.add<ValueArg<int>>(opt.verbosity, "", "verbosity",
        "verbosity level (0:silent,1:quiet,2:improvements only,3:verbose",
        false, 2, "int");

    cmd.add<ValueArg<int>>(opt.randwalkiter, "", "randwalkiter",
        "number of random walk iterations during local search", false, 0,
        "int");
    cmd.add<ValueArg<long int>>(opt.lsiter, "", "lsiter",
        "number of local search iterations", false, 0, "int");
    cmd.add<ValueArg<int>>(opt.lsextra, "", "lsextra",
        "number of local search extra iterations during the i-dsatur phase",
        false, 0, "int");
    cmd.add<ValueArg<int>>(opt.dsatlimit, "", "dsatlimit",
        "iteration limit during dsat moves", false, 0, "int");
    cmd.add<SwitchArg>(opt.switchdescent, "", "switchdescent",
        "use descent move in LS", false);
    cmd.add<SwitchArg>(opt.switchreact, "", "switchreact",
        "do (not) use react move in LS", true);
    cmd.add<SwitchArg>(
        opt.focus, "", "focus", "do not move the core in LS", false);
    cmd.add<ValueArg<int>>(
        opt.rw, "", "rw", "frequency of path exploration", false, 1, "int");

    cmd.add<ValueArg<int>>(
        opt.seed, "", "seed", "random seed", false, 12345, "int");
    cmd.add<ValueArg<int>>(
        opt.tenure, "", "tenure", "tabu tenure", false, 10, "int");

    cmd.add<ValueArg<int>>(opt.dynamiclimit, "", "dynamiclimit",
        "update the LS iteration limit dynamically [0=static, 1=geometric, 2=backoff]", false, 2, "int");
    cmd.add<SwitchArg>(opt.dynrandpath, "", "dynrandpath",
        "update the randpath ratio dynamically", true);

    cmd.add<ValueArg<int>>(opt.dynfactor, "", "dynfactor",
        "geometric factor for iter limit", false, 2, "int");
    cmd.add<ValueArg<int>>(opt.dyndiv, "", "dyndiv",
        "geometric divisor for iter limit", false, 1, "int");

    cmd.add<ValueArg<int>>(opt.rpfactor, "", "rpfactor",
        "factor for randpath ratio update", false, 2, "int");
    cmd.add<ValueArg<int>>(opt.rpdiv, "", "rpdiv",
        "divisor for randpath ratio update", false, 1, "int");
    cmd.add<ValueArg<int>>(opt.rpmin, "", "rpmin",
        "min for randpath ratio update", false, 3, "int");
    cmd.add<ValueArg<int>>(opt.rpmax, "", "rpmax",
        "max for randpath ratio update", false, 20, "int");

    cmd.add<SwitchArg>(opt.norecolor, "", "norecolor",
        "do not use recolor in dsatur", false);

    cmd.parse(argc, argv);
    return opt;
}

void options::describe(std::ostream& os)
{
    os << "[options] GC configuration\n";
    os << "[options] cmdline = " << cmdline << "\n";
    os << "[options] Instance file = " << instance_file << "\n";
    os << "[options] Clause learning = " << learning << "\n";
    os << "[options] Polarity policy = " << polarity << "\n";
    os << "[options] Clique ordering = " << ordering << "\n";
    os << "[options]  ... low degree = " << ordering_low_degree << "\n";
    os << "[options] Color variables = " << (xvars ? "present" : "absent") << "\n";
    os << "[options] Bound policy    = " << boundalg << "\n";
    os << "[options] Adaptive bounds = " << adaptive << "\n";
    os << "[options] Preprocessing   = " << (preprocessing == 0 ? "none" : preprocessing == 1 ?  "dense" : "sparse") << "\n";
    os << "[options] IS constraints  = " << (int)indset_constraints << "\n";
    os << "[options] fillin          = "
       << (fillin == options::FILLIN_NONE
                  ? "no"
                  : (fillin == options::FILLIN_PROPAGATE ? "yes (propagate)"
                                                         : "yes (decompose)"))
       << "\n";
    os << "[options] Triangle UP     = " << triangle_up << "\n";
    os << "[options] Branching strat = ";
    switch (branching) {
    case gc::options::VSIDS:
        os << "VSIDS\n";
        break;
    case gc::options::BRELAZ:
        os << "Brelaz\n";
        break;
    case gc::options::PARTITION_PRODUCT:
        os << "Largest edge partition size (x)\n";
        break;
    case gc::options::PARTITION_SUM:
        os << "Largest edge partition size (+)\n";
        break;
    case gc::options::DEGREE_PRODUCT:
        os << "Largest edge degree (x)\n";
        break;
    case gc::options::DEGREE_SUM:
        os << "Largest edge degree (+)\n";
        break;
    case gc::options::DEGREE_UNION:
        os << "Largest edge neighborhood\n";
        break;
    case gc::options::PARTITION_PRODUCT_DYN:
        os << "Largest edge partition size (x) (select)\n";
        break;
    case gc::options::PARTITION_SUM_DYN:
        os << "Largest edge partition size (+) (select)\n";
        break;
    case gc::options::DEGREE_PRODUCT_DYN:
        os << "Largest edge degree (x) (select)\n";
        break;
    case gc::options::DEGREE_SUM_DYN:
        os << "Largest edge degree (+) (select)\n";
        break;
    case gc::options::DEGREE_UNION_DYN:
        os << "Largest edge neighborhood (select)\n";
        break;
    case gc::options::VSIDS_PHASED:
        os << "VSIDS with partition phase heuristic\n";
        break;
    case gc::options::VSIDS_GUIDED:
        os << "VSIDS with solution phase saving\n";
        break;
    case gc::options::VSIDS_CLIQUE:
        os << "VSIDS restricted to merging with large cliques\n";
        break;
    case gc::options::VSIDS_COLORS_POSITIVE:
        os << "VSIDS restricted to assignments to color variables\n";
        break;
    case gc::options::VERTEX_ACTIVITY:
        os << "Exponentially decaying vertex activity\n";
        break;
    case gc::options::VERTEX_DOM_OVER_ACT:
        os << "Domain size over Exponentially decaying vertex activity\n";
        break;
    }
    os << "[options]  ... low degree = " << branching_low_degree << "\n";
    os << "[options]  ... brelaz fst = " << brelaz_first << "\n";
    os << "[options]  ... phase save = " << phase_saving << "\n";
    os << "[options] Strategy        = ";
    switch (strategy) {
    case BNB:
        os << "branch and bound";
        break;
    case BOTTOMUP:
        os << "bottom-up";
        break;
    case LOCALSEARCH:
        os << "local search";
        break;
    case BOUNDS:
        os << "bounds";
        break;
    case CLEVER:
        os << "(Verma et al.)";
        break;
    case COLOR:
        os << "heuristic coloring";
        break;
    case TEST:
        os << "test";
        break;
    case IDSATUR:
        os << "iterated dsatur";
        break;
    }
    os << "\n";
    os << std::endl;
}

} // namespace gc
