set_languages("c++17")

includes("../common/xmake.lua")


target("reservoir")
    set_kind("static")
    add_headerfiles("*.h")
    add_files("*.cpp")
    add_includedirs(".",{public=true})
    add_deps("common")
