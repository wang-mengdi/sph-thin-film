add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_languages("c++17")

add_requires("perlinnoise >=3.0.0")
add_requires("libomp")
add_requires("tinyobjloader v2.0.0rc13")


set_rundir("$(projectdir)")

includes("./src/physics/xmake.lua")

target("sph_bubble")
    set_kind("binary")
    add_headerfiles("sph_bubble/*.h")
    add_files("sph_bubble/*.cpp")
    add_includedirs("sph_bubble", {public = true})
    add_packages("perlinnoise", {public = true})
    add_packages("libomp")
    add_ldflags("-fopenmp")
    add_deps("physics")