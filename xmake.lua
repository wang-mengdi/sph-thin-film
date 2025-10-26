add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_languages("c++17")

if is_plat("windows") then
    add_defines("WIN32")
elseif is_plat("linux") then
    add_defines("__linux__")
end


add_requires("perlinnoise >=3.0.0")
add_requires("libomp")
add_requires("tinyobjloader v2.0.0rc13")
add_requires("freeglut")
add_requires("glew")
add_requires("glm")
add_requires("stb")
add_requires("imgui")


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

target("viewer")
    set_kind("binary")
    add_headerfiles("viewer/*.h")
    add_files("viewer/*.cpp")
    add_includedirs("viewer", {public = true})
    add_packages("freeglut")
    add_packages("glew")
    add_packages("stb")
    add_packages("glm")
    add_packages("imgui")
    add_deps("geometry")    