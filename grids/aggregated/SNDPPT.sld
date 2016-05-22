<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>SNDPPT</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
          <Geometry>
              <PropertyName>Sand content (50-2000 micro meter) mass fraction in %</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="intervals">
              <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
              <ColorMapEntry color="#FFFF00" quantity="1" label="0 - 1" opacity="0.7"/>
              <ColorMapEntry color="#F8F806" quantity="3" label="1 - 3" opacity="0.7"/>
              <ColorMapEntry color="#F1F10C" quantity="4" label="3 - 4" opacity="0.7"/>
              <ColorMapEntry color="#EBEB13" quantity="6" label="4 - 6" opacity="0.7"/>
              <ColorMapEntry color="#E4E419" quantity="8" label="6 - 8" opacity="0.7"/>
              <ColorMapEntry color="#DDDD20" quantity="10" label="8 - 10" opacity="0.7"/>
              <ColorMapEntry color="#D7D726" quantity="12" label="10 - 12" opacity="0.7"/>
              <ColorMapEntry color="#D0D02D" quantity="14" label="12 - 14" opacity="0.7"/>
              <ColorMapEntry color="#CACA33" quantity="16" label="14 - 16" opacity="0.7"/>
              <ColorMapEntry color="#C3C33A" quantity="19" label="16 - 19" opacity="0.7"/>
              <ColorMapEntry color="#BCBC41" quantity="21" label="19 - 21" opacity="0.7"/>
              <ColorMapEntry color="#B6B647" quantity="23" label="21 - 23" opacity="0.7"/>
              <ColorMapEntry color="#B0B04E" quantity="26" label="23 - 26" opacity="0.7"/>
              <ColorMapEntry color="#A9A954" quantity="28" label="26 - 28" opacity="0.7"/>
              <ColorMapEntry color="#A3A35A" quantity="30" label="28 - 30" opacity="0.7"/>
              <ColorMapEntry color="#9C9C61" quantity="32" label="30 - 32" opacity="0.7"/>
              <ColorMapEntry color="#959568" quantity="35" label="32 - 35" opacity="0.7"/>
              <ColorMapEntry color="#8F8F6E" quantity="37" label="35 - 37" opacity="0.7"/>
              <ColorMapEntry color="#898975" quantity="39" label="37 - 39" opacity="0.7"/>
              <ColorMapEntry color="#82827B" quantity="41" label="39 - 41" opacity="0.7"/>
              <ColorMapEntry color="#7B7B82" quantity="44" label="41 - 44" opacity="0.7"/>
              <ColorMapEntry color="#757589" quantity="46" label="44 - 46" opacity="0.7"/>
              <ColorMapEntry color="#6E6E8F" quantity="48" label="46 - 48" opacity="0.7"/>
              <ColorMapEntry color="#686895" quantity="51" label="48 - 51" opacity="0.7"/>
              <ColorMapEntry color="#61619C" quantity="53" label="51 - 53" opacity="0.7"/>
              <ColorMapEntry color="#5A5AA3" quantity="56" label="53 - 56" opacity="0.7"/>
              <ColorMapEntry color="#5454A9" quantity="58" label="56 - 58" opacity="0.7"/>
              <ColorMapEntry color="#4D4DB0" quantity="61" label="58 - 60.4" opacity="0.7"/>
              <ColorMapEntry color="#4747B6" quantity="63" label="60.4 - 63" opacity="0.7"/>
              <ColorMapEntry color="#4141BC" quantity="66" label="63 - 66" opacity="0.7"/>
              <ColorMapEntry color="#3A3AC3" quantity="68" label="66 - 68" opacity="0.7"/>
              <ColorMapEntry color="#3333CA" quantity="71" label="68 - 71" opacity="0.7"/>
              <ColorMapEntry color="#2D2DD0" quantity="74" label="71 - 74" opacity="0.7"/>
              <ColorMapEntry color="#2626D7" quantity="77" label="74 - 77" opacity="0.7"/>
              <ColorMapEntry color="#2020DD" quantity="80" label="77 - 80" opacity="0.7"/>
              <ColorMapEntry color="#1919E4" quantity="83" label="80 - 83" opacity="0.7"/>
              <ColorMapEntry color="#1212EB" quantity="86" label="83 - 86" opacity="0.7"/>
              <ColorMapEntry color="#0C0CF1" quantity="90" label="86 - 90" opacity="0.7"/>
              <ColorMapEntry color="#0606F8" quantity="94" label="90 - 94" opacity="0.7"/>
              <ColorMapEntry color="#0000FF" quantity="100" label="94 - 100" opacity="0.7"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>