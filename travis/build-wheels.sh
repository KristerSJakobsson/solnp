#!/bin/bash
set -e -x

yum install -y cmake
yum install -y openssl-devel

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
    "${PYBIN}/python" /io/setup.py sdist -d /io/wheelhouse/
done

# Bundle external shared libraries into the wheels
for WHEEL in wheelhouse/*.whl; do
    auditwheel repair "$WHEEL" -w /io/wheelhouse/
done

# Test
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/requirements-dev.txt
    "${PYBIN}/pip" install --no-index -f /io/wheelhouse pysolnp
    (cd "$PYHOME"; "${PYBIN}/pytest" /io/python_solnp/test/test.py)
done

#  Upload
/opt/python/cp37-cp37m/bin/pip install twine==2.3.2  # Note: Old version to provide support with i686 anylinux docker
for WHEEL in /io/wheelhouse/pysolnp*; do
    /opt/python/cp37-cp37m/bin/twine upload \
        --skip-existing \
        --repository-url "${TWINE_SERVER}" \
        --username "${TWINE_USERNAME}" \
        --password "${TWINE_PASSWORD}" \
        "${WHEEL}"
done