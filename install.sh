python setup.py sdist
cd dist
tar zxvf *.tar.gz
cd PyGS*
sudo python setup.py install
cd ..
sudo rm -rf PyGS*/*
rmdir PyGS*
cd ..