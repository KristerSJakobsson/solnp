docker run --rm -v $(pwd):/io -e TWINE_USERNAME="$PYPI_USERNAME" -e TWINE_PASSWORD="$PYPI_PASSWORD" -e TWINE_SERVER="$PYPI_SERVER" "$DOCKER_IMAGE" "$PRE_CMD" /io/travis/build-wheels.sh
ls wheelhouse/
