include_directories("../../src")

add_cnoid_simple_controller(my_controller  main.cpp myrobot.cpp myrobot.h)

target_link_libraries(my_controller vnoid_lib)
