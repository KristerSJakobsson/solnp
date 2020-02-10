#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y atlas-devel
yum install -y cmake

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/requirements-dev.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
    "${PYBIN}/python" /io/setup.py sdist -d /io/wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

for whl in /io/wheelhouse/*; do
    echo "$whl"
done

# Install package and test
for PYBIN in /opt/python/*/bin; do
#    ln -s "${PYBIN}/cmake" /usr/bin/cmake
    "${PYBIN}/pip" install pysolnp --no-index -f /io/wheelhouse
    (cd "$PYHOME"; "${PYBIN}/nosetests" /io/python_solnp/test/test.py)
done
#
##  Upload
#for WHEEL in /io/wheelhouse/pysolnp*; do
#    /opt/python/cp37-cp37m/bin/twine upload \
#        --skip-existing \
#        --repository-url "${TWINE_SERVER}" \
#        --username "${TWINE_USERNAME}" \
#        --password "${TWINE_PASSWORD}" \
#        "${WHEEL}"
#done