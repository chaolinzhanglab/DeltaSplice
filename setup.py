from setuptools import setup
import io


setup(name='deltasplice',
      description='DeltaSplice: A neural network model to predict splice site usage and splicing-altering mutations',
      long_description=io.open('README.md', encoding='utf-8').read(),
      long_description_content_type='text/markdown',
      version='1.0.0',
      packages=['deltasplice', "deltasplice.models"],
      package_data={'deltasplice': [
            'pretrained_models/DeltaSplice_models/model.ckpt-0',
            'pretrained_models/DeltaSplice_models/model.ckpt-1',
            'pretrained_models/DeltaSplice_models/model.ckpt-2',
            'pretrained_models/DeltaSplice_models/model.ckpt-3',
            'pretrained_models/DeltaSplice_models/model.ckpt-4',
            "annotations/grch37.txt",
            "annotations/grch38.txt",
            "annotations/",
            "data/anno/data.json"
                                ]},
      entry_points={'console_scripts': ['deltasplice=deltasplice.__main__:main']})