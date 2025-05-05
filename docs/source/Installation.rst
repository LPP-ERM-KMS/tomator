Installation
============

Cloning & Building
------------------

First, clone the repository to a location you'd like::

    git clone https://github.com/LPP-ERM-KMS/tomator.git

Then navigate to the src folder in tomator::

    cd tomator/src

and build the software (replace the number after j with the number you get when running nproc)::

    make -j8

This will have built an executable called 'Tomator1D' which is the primary
binary, you are now done and may move on to :doc:`Usage`.  If modifications were made
to the software or you wish to rebuild for other reasons make sure to clean
first before rebuilding::

    make clean && make -j8


Troubleshooting
---------------

arm-macs are not yet supported (see https://github.com/LPP-ERM-KMS/tomator/issues/5)
