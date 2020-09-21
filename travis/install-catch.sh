# Custom shell script for installing catch on travis because the apt install path doesn't work with CMake
sudo apt-get -y install catch
ln -s /usr/include/catch/catch.hpp /usr/include/catch.hpp
