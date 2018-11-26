from setuptools import setup, find_packages
import versioneer

setup(name='SearchSifter',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          "pymysql",
      ],
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      )
