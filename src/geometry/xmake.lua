set_languages("c++17")

includes("../reservoir/xmake.lua")

add_requires("nanoflann =1.3.2", {configs = {cxx11 = true}})

target("geometry")
    set_kind("static")
    add_headerfiles("*.h")
    add_files("*.cpp")
    add_includedirs(".",{public=true})
    add_packages("nanoflann", {public = true})
    add_deps("reservoir")
