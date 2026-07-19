#include "xpclr.hpp"

int main(int argc, char** argv) {
    xpclr::Options opt = xpclr::parse_args(argc, argv);
    return xpclr::run_xpclr(opt);
}
