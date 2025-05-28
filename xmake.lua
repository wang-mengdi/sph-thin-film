add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_languages("c++17")

add_requires("perlinnoise >=3.0.0")
add_requires("libomp >=19.1.0")
add_requires("tinyobjloader v2.0.0rc13")


set_rundir("$(projectdir)")

includes("./src/geometry/xmake.lua")

target("melp")
    set_kind("binary")
    add_headerfiles("melp/*.h")
    add_files("melp/*.cpp")
    add_includedirs("melp", {public = true})
    add_packages("perlinnoise", {public = true})
    add_packages("libomp")
    add_ldflags("-fopenmp")
    add_deps("geometry")