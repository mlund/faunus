from conans import ConanFile, CMake

class FaunusConan(ConanFile):
    generators = "cmake"
    settings = "os", "compiler", "build_type"
    requires = "cereal/1.3.0", "docopt.cpp/0.6.3", "doctest/2.4.4", "eigen/3.3.9", "exprtk/20181202@bincrafters/stable", "fmt/7.1.3", "nanobench/4.3.0", "nlohmann_json/3.9.1", "pybind11/2.6.1", "range-v3/0.11.0", "spdlog/1.8.2", "tbb/2020.3", "trompeloeil/39"
    default_options = {"nlohmann_json:multiple_headers": True, "tbb:shared": False}

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build(target="faunus")
        cmake.build(target="pyfaunus")
        cmake.test()
