<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>ORCDRC</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
           <Geometry>
              <PropertyName>Soil organic carbon content (fine earth fraction) in g per kg</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" quantity="-32768" label="NODATA" opacity="0.0"/>
              <ColorMapEntry color="#000180" quantity="0.2" label="0 - 0.2" opacity="0.7"/>
              <ColorMapEntry color="#000393" quantity="0.4" label="0.2 - 0.4" opacity="0.7"/>
              <ColorMapEntry color="#0006A6" quantity="0.6" label="0.4 - 0.6" opacity="0.7"/>
              <ColorMapEntry color="#000FB7" quantity="0.8" label="0.6 - 0.8" opacity="0.7"/>
              <ColorMapEntry color="#0018C8" quantity="1.1" label="0.8 - 1.1" opacity="0.7"/>
              <ColorMapEntry color="#0025D6" quantity="1.5" label="1.1 - 1.5" opacity="0.7"/>
              <ColorMapEntry color="#0032E3" quantity="1.9" label="1.5 - 1.9" opacity="0.7"/>
              <ColorMapEntry color="#0043ED" quantity="2.4" label="1.9 - 2.4" opacity="0.7"/>
              <ColorMapEntry color="#0054F6" quantity="3.0" label="2.4 - 3" opacity="0.7"/>
              <ColorMapEntry color="#0067FA" quantity="3.6" label="3 - 3.6" opacity="0.7"/>
              <ColorMapEntry color="#0079FE" quantity="4.4" label="3.6 - 4.4" opacity="0.7"/>
              <ColorMapEntry color="#038DFC" quantity="5.3" label="4.4 - 5.3" opacity="0.7"/>
              <ColorMapEntry color="#06A0F9" quantity="6.3" label="5.3 - 6.3" opacity="0.7"/>
              <ColorMapEntry color="#0DB2F2" quantity="7.5" label="6.3 - 7.5" opacity="0.7"/>
              <ColorMapEntry color="#15C3E9" quantity="8.9" label="7.5 - 8.9" opacity="0.7"/>
              <ColorMapEntry color="#21D2DD" quantity="10.5" label="8.9 - 10.5" opacity="0.7"/>
              <ColorMapEntry color="#2FE0CF" quantity="12.4" label="10.5 - 12.4" opacity="0.7"/>
              <ColorMapEntry color="#3EEBC0" quantity="14.6" label="12.4 - 14.6" opacity="0.7"/>
              <ColorMapEntry color="#4FF3AF" quantity="17.2" label="14.6 - 17.2" opacity="0.7"/>
              <ColorMapEntry color="#62F99C" quantity="20.2" label="17.2 - 20.2" opacity="0.7"/>
              <ColorMapEntry color="#75FD8A" quantity="23.7" label="20.2 - 23.7" opacity="0.7"/>
              <ColorMapEntry color="#89FD76" quantity="27.8" label="23.7 - 27.8" opacity="0.7"/>
              <ColorMapEntry color="#9BF963" quantity="32.6" label="27.8 - 32.6" opacity="0.7"/>
              <ColorMapEntry color="#AEF450" quantity="38.1" label="32.6 - 38.1" opacity="0.7"/>
              <ColorMapEntry color="#BFEB3F" quantity="44.6" label="38.1 - 44.6" opacity="0.7"/>
              <ColorMapEntry color="#CFE02F" quantity="52.1" label="44.6 - 52.1" opacity="0.7"/>
              <ColorMapEntry color="#DCD322" quantity="60.9" label="52.1 - 60.9" opacity="0.7"/>
              <ColorMapEntry color="#E8C416" quantity="71.1" label="60.9 - 71.1" opacity="0.7"/>
              <ColorMapEntry color="#F1B30D" quantity="83.0" label="71.1 - 83" opacity="0.7"/>
              <ColorMapEntry color="#F8A106" quantity="96.9" label="83 - 96.9" opacity="0.7"/>
              <ColorMapEntry color="#FC8F03" quantity="113.0" label="96.9 - 113" opacity="0.7"/>
              <ColorMapEntry color="#FE7B00" quantity="131.8" label="113 - 132" opacity="0.7"/>
              <ColorMapEntry color="#FA6800" quantity="153.8" label="132 - 154" opacity="0.7"/>
              <ColorMapEntry color="#F65500" quantity="179.3" label="154 - 179" opacity="0.7"/>
              <ColorMapEntry color="#ED4400" quantity="209.1" label="179 - 209" opacity="0.7"/>
              <ColorMapEntry color="#E43300" quantity="243.8" label="209 - 244" opacity="0.7"/>
              <ColorMapEntry color="#D62500" quantity="284.2" label="244 - 284" opacity="0.7"/>
              <ColorMapEntry color="#C91800" quantity="331.3" label="284 - 331" opacity="0.7"/>
              <ColorMapEntry color="#B80F00" quantity="386.1" label="331 - 386" opacity="0.7"/>
              <ColorMapEntry color="#A70700" quantity="450.0" label="386 - 450" opacity="0.7"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>