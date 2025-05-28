set_languages("c++17")

includes("../geometry/xmake.lua")

target("physics")
    set_kind("static")
    add_headerfiles("*.h")
    add_files("*.cpp")
    add_includedirs(".",{public=true})
    add_deps("geometry")
