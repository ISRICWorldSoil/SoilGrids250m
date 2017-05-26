<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>SoilGrids250m_TAXOUSDA</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
                     <Geometry>
              <ogc:PropertyName>USDA 2014 class</ogc:PropertyName>
            </Geometry>
            <ColorMap type="values">
              <ColorMapEntry color="#00FF00" quantity="70" label="Albolls" opacity="1.0"/>
              <ColorMapEntry color="#ADFF2F" quantity="80" label="Aqualfs" opacity="1.0"/>
              <ColorMapEntry color="#FF00FF" quantity="20" label="Aquands" opacity="1.0"/>
              <ColorMapEntry color="#7FFFD4" quantity="95" label="Aquents" opacity="1.0"/>
              <ColorMapEntry color="#CF595C" quantity="90" label="Aquepts" opacity="1.0"/>
              <ColorMapEntry color="#FFFF00" quantity="40" label="Aquerts" opacity="1.0"/>
              <ColorMapEntry color="#D8BFD8" quantity="15" label="Aquods" opacity="1.0"/>
              <ColorMapEntry color="#03FF05" quantity="71" label="Aquolls" opacity="1.0"/>
              <ColorMapEntry color="#FF0000" quantity="30" label="Aquox" opacity="1.0"/>
              <ColorMapEntry color="#FFA500" quantity="60" label="Aquults" opacity="1.0"/>
              <ColorMapEntry color="#7DFFD2" quantity="96" label="Arents" opacity="1.0"/>
              <ColorMapEntry color="#FFDDC2" quantity="54" label="Argids" opacity="1.0"/>
              <ColorMapEntry color="#09FE03" quantity="69" label="Borolls" opacity="1.0"/>
              <ColorMapEntry color="#E7CDC0" quantity="55" label="Calcids" opacity="1.0"/>
              <ColorMapEntry color="#F3E3C8" quantity="56" label="Cambids" opacity="1.0"/>
              <ColorMapEntry color="#A5FF2F" quantity="81" label="Cryalfs" opacity="1.0"/>
              <ColorMapEntry color="#FA02FA" quantity="21" label="Cryands" opacity="1.0"/>
              <ColorMapEntry color="#E05C5D" quantity="92" label="Cryepts" opacity="1.0"/>
              <ColorMapEntry color="#FFDAB9" quantity="50" label="Cryids" opacity="1.0"/>
              <ColorMapEntry color="#D4C3D4" quantity="16" label="Cryods" opacity="1.0"/>
              <ColorMapEntry color="#0FEA03" quantity="74" label="Cryolls" opacity="1.0"/>
              <ColorMapEntry color="#F5D3B9" quantity="52" label="Durids" opacity="1.0"/>
              <ColorMapEntry color="#B22328" quantity="11" label="Fibrists" opacity="1.0"/>
              <ColorMapEntry color="#73FFD2" quantity="98" label="Fluvents" opacity="1.0"/>
              <ColorMapEntry color="#A52A2A" quantity="10" label="Folists" opacity="1.0"/>
              <ColorMapEntry color="#EB05EB" quantity="27" label="Gelands" opacity="1.0"/>
              <ColorMapEntry color="#CB5A5F" quantity="86" label="Gelepts" opacity="1.0"/>
              <ColorMapEntry color="#DDB9DD" quantity="19" label="Gelods" opacity="1.0"/>
              <ColorMapEntry color="#E8C8B8" quantity="53" label="Gypsids" opacity="1.0"/>
              <ColorMapEntry color="#B41919" quantity="12" label="Hemists" opacity="1.0"/>
              <ColorMapEntry color="#48D1CC" quantity="5" label="Histels" opacity="1.0"/>
              <ColorMapEntry color="#D2BAD2" quantity="17" label="Humods" opacity="1.0"/>
              <ColorMapEntry color="#F3A702" quantity="61" label="Humults" opacity="1.0"/>
              <ColorMapEntry color="#CA5960" quantity="89" label="Ochrepts" opacity="1.0"/>
              <ColorMapEntry color="#4EC8CC" quantity="7" label="Orthels" opacity="1.0"/>
              <ColorMapEntry color="#88EEC8" quantity="99" label="Orthents" opacity="1.0"/>
              <ColorMapEntry color="#D5C0D5" quantity="18" label="Orthods" opacity="1.0"/>
              <ColorMapEntry color="#FB0202" quantity="33" label="Perox" opacity="1.0"/>
              <ColorMapEntry color="#86F5CD" quantity="97" label="Psamments" opacity="1.0"/>
              <ColorMapEntry color="#05F300" quantity="72" label="Rendolls" opacity="1.0"/>
              <ColorMapEntry color="#F5D7BB" quantity="51" label="Salids" opacity="1.0"/>
              <ColorMapEntry color="#A42828" quantity="13" label="Saprists" opacity="1.0"/>
              <ColorMapEntry color="#FC05FA" quantity="22" label="Torrands" opacity="1.0"/>
              <ColorMapEntry color="#EBEB0C" quantity="43" label="Torrerts" opacity="1.0"/>
              <ColorMapEntry color="#F50505" quantity="31" label="Torrox" opacity="1.0"/>
              <ColorMapEntry color="#43D4D2" quantity="6" label="Turbels" opacity="1.0"/>
              <ColorMapEntry color="#8CFF19" quantity="84" label="Udalfs" opacity="1.0"/>
              <ColorMapEntry color="#F100F1" quantity="26" label="Udands" opacity="1.0"/>
              <ColorMapEntry color="#CD5C5C" quantity="85" label="Udepts" opacity="1.0"/>
              <ColorMapEntry color="#EEFF06" quantity="45" label="Uderts" opacity="1.0"/>
              <ColorMapEntry color="#0CFF0C" quantity="76" label="Udolls" opacity="1.0"/>
              <ColorMapEntry color="#FF0E0E" quantity="34" label="Udox" opacity="1.0"/>
              <ColorMapEntry color="#FB9C00" quantity="62" label="Udults" opacity="1.0"/>
              <ColorMapEntry color="#8CFF37" quantity="82" label="Ustalfs" opacity="1.0"/>
              <ColorMapEntry color="#F50CF0" quantity="25" label="Ustands" opacity="1.0"/>
              <ColorMapEntry color="#D35740" quantity="93" label="Ustepts" opacity="1.0"/>
              <ColorMapEntry color="#F5EB00" quantity="44" label="Usterts" opacity="1.0"/>
              <ColorMapEntry color="#00F000" quantity="75" label="Ustolls" opacity="1.0"/>
              <ColorMapEntry color="#F20A0A" quantity="32" label="Ustox" opacity="1.0"/>
              <ColorMapEntry color="#F0B005" quantity="63" label="Ustults" opacity="1.0"/>
              <ColorMapEntry color="#FC04F5" quantity="24" label="Vitrands" opacity="1.0"/>
              <ColorMapEntry color="#AFFF19" quantity="83" label="Xeralfs" opacity="1.0"/>
              <ColorMapEntry color="#AFFF08" quantity="23" label="Xerands" opacity="1.0"/>
              <ColorMapEntry color="#D95F35" quantity="94" label="Xerepts" opacity="1.0"/>
              <ColorMapEntry color="#FAFA05" quantity="66" label="Xererts" opacity="1.0"/>
              <ColorMapEntry color="#02F00A" quantity="73" label="Xerolls" opacity="1.0"/>
              <ColorMapEntry color="#F7980F" quantity="64" label="Xerults" opacity="1.0"/>
              <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>
