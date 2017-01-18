#!groovy
node {

    stage 'checkout'
    git url: 'https://github.com/justincely/lightcurve'

    stage 'setup'
    sh 'wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;'
    sh 'chmod +x ./miniconda.sh'
    sh 'ls -la'

    sh 'echo "home"'
    sh 'echo $HOME'
    sh 'echo "path"'
    sh 'echo $PATH'

    sh './miniconda.sh -b -p $HOME/miniconda'
    sh 'export PATH="$HOME/miniconda/bin:$PATH"'

    sh 'echo $PATH'

    // list conda environments
    sh 'echo "checking conda envs"'
    sh 'ls $HOME/miniconda/bin'
    sh 'conda env list'

    stage 'install'
    sh 'python setup.py install'

    stage 'build'
    sh 'nosetests'


}
