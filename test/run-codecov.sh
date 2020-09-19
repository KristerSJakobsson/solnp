cmake .
make gcov
make lcov
bash <(curl -s https://codecov.io/bash) -X gcov || echo "Codecov did not collect coverage reports"
