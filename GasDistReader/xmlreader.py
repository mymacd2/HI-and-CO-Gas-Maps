import xml.etree.ElementTree as ET

tree = ET.parse('HI_model.xml')
root = tree.getroot()


for element in root.iter('RadialProfile'):
    print(element.attrib)

