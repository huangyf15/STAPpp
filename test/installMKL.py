import platform
import os

assert platform.system() == 'Linux'

if os.path.exists('/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/include') and \
   os.path.exists('/opt/intel/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64') :
    print('MKL installed or cached.')
else:
    print('installing MKL.')
    f = os.path.dirname
    projectDir = f(f(os.path.abspath(__file__)))
    os.chdir(projectDir + os.sep + 'temp')
    os.system(
        'wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB')
    os.system('sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB')
    os.system(
        'sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list')
    os.system('sudo apt-get update')
    os.system('sudo apt-get install intel-mkl-2018.1-038 -y')
    print('MKL installed.')
