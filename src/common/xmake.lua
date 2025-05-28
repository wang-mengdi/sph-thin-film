set_languages("c++17")

add_requires("fmt =8.1.1")
add_requires("nlohmann_json >=3.12.0")
add_requires("eigen >=3.4.0")


target("common")
    set_kind("static")
    add_headerfiles("*.h", "physbam/*.h")
    add_files("*.cpp")
    add_includedirs(".","physbam",{public=true})
    add_packages("fmt",{public=true})
    add_packages("nlohmann_json",{public=true})
    add_packages("eigen",{public=true})
