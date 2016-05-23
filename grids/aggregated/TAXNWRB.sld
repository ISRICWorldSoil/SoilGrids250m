<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
  <NamedLayer>
    <Name>SoilGrids250m</Name>
    <UserStyle>
      <Title>TAXNWRB</Title>
      <FeatureTypeStyle>
        <Rule>
          <RasterSymbolizer>
           <Geometry>
              <PropertyName>WRB 2006 class</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
            <ColorMap type="values">
              <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
              <ColorMapEntry color="#FC894D" quantity="43" label="Acric Ferralsols" opacity="0.7"/>
              <ColorMapEntry color="#845851" quantity="95" label="Acric Plinthosols" opacity="0.7"/>
              <ColorMapEntry color="#FFE8BE" quantity="15" label="Albic Arenosols" opacity="0.7"/>
              <ColorMapEntry color="#F896B4" quantity="76" label="Albic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#A05246" quantity="85" label="Alic Nitisols" opacity="0.7"/>
              <ColorMapEntry color="#FC6B5D" quantity="12" label="Aluandic Andosols" opacity="0.7"/>
              <ColorMapEntry color="#FDE4B7" quantity="99" label="Aric Regosols" opacity="0.7"/>
              <ColorMapEntry color="#FDD2B8" quantity="100" label="Calcaric Regosols" opacity="0.7"/>
              <ColorMapEntry color="#E5C75D" quantity="36" label="Calcic Chernozems" opacity="0.7"/>
              <ColorMapEntry color="#AB9ABF" quantity="53" label="Calcic Gleysols" opacity="0.7"/>
              <ColorMapEntry color="#FFF9AD" quantity="59" label="Calcic Gypsisols" opacity="0.7"/>
              <ColorMapEntry color="#8291B8" quantity="61" label="Calcic Histosols" opacity="0.7"/>
              <ColorMapEntry color="#D2A18D" quantity="66" label="Calcic Kastanozems" opacity="0.7"/>
              <ColorMapEntry color="#F9A890" quantity="77" label="Calcic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#F8D8ED" quantity="108" label="Calcic Solonetz" opacity="0.7"/>
              <ColorMapEntry color="#A87188" quantity="115" label="Calcic Vertisols" opacity="0.7"/>
              <ColorMapEntry color="#8A9698" quantity="62" label="Cryic Histosols" opacity="0.7"/>
              <ColorMapEntry color="#E4DBBF" quantity="10" label="Cutanic Alisols" opacity="0.7"/>
              <ColorMapEntry color="#FDCE70" quantity="25" label="Endogleyic Cambisols" opacity="0.7"/>
              <ColorMapEntry color="#D7965F" quantity="90" label="Endogleyic Planosols" opacity="0.7"/>
              <ColorMapEntry color="#FEE7C0" quantity="16" label="Ferralic Arenosols" opacity="0.7"/>
              <ColorMapEntry color="#EDC669" quantity="26" label="Ferralic Cambisols" opacity="0.7"/>
              <ColorMapEntry color="#788092" quantity="63" label="Fibric Histosols" opacity="0.7"/>
              <ColorMapEntry color="#F9A091" quantity="78" label="Gleyic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#6EAA7C" quantity="97" label="Gleyic Podzols" opacity="0.7"/>
              <ColorMapEntry color="#E2CDCC" quantity="109" label="Gleyic Solonetz" opacity="0.7"/>
              <ColorMapEntry color="#E05B9A" quantity="105" label="Gypsic Solonchaks" opacity="0.7"/>
              <ColorMapEntry color="#FE813E" quantity="1" label="Haplic Acrisols" opacity="0.7"/>
              <ColorMapEntry color="#FD9F39" quantity="2" label="Haplic Acrisols (Alumic)" opacity="0.7"/>
              <ColorMapEntry color="#FDAE6B" quantity="3" label="Haplic Acrisols (Ferric)" opacity="0.7"/>
              <ColorMapEntry color="#FD8D3C" quantity="4" label="Haplic Acrisols (Humic)" opacity="0.7"/>
              <ColorMapEntry color="#F0C85C" quantity="7" label="Haplic Albeluvisols" opacity="0.7"/>
              <ColorMapEntry color="#F5EBCC" quantity="11" label="Haplic Alisols" opacity="0.7"/>
              <ColorMapEntry color="#F5E7CA" quantity="13" label="Haplic Andosols" opacity="0.7"/>
              <ColorMapEntry color="#FEE3C0" quantity="17" label="Haplic Arenosols" opacity="0.7"/>
              <ColorMapEntry color="#FEE4B1" quantity="18" label="Haplic Arenosols (Calcaric)" opacity="0.7"/>
              <ColorMapEntry color="#FFEF51" quantity="21" label="Haplic Calcisols" opacity="0.7"/>
              <ColorMapEntry color="#F8E729" quantity="22" label="Haplic Calcisols (Sodic)" opacity="0.7"/>
              <ColorMapEntry color="#FDE260" quantity="27" label="Haplic Cambisols" opacity="0.7"/>
              <ColorMapEntry color="#FDE770" quantity="28" label="Haplic Cambisols (Calcaric)" opacity="0.7"/>
              <ColorMapEntry color="#FDF770" quantity="29" label="Haplic Cambisols (Chromic)" opacity="0.7"/>
              <ColorMapEntry color="#F6CE4B" quantity="30" label="Haplic Cambisols (Dystric)" opacity="0.7"/>
              <ColorMapEntry color="#FDE170" quantity="31" label="Haplic Cambisols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#FDE05F" quantity="32" label="Haplic Cambisols (Humic)" opacity="0.7"/>
              <ColorMapEntry color="#FDEB70" quantity="33" label="Haplic Cambisols (Sodic)" opacity="0.7"/>
              <ColorMapEntry color="#E5D15C" quantity="37" label="Haplic Chernozems" opacity="0.7"/>
              <ColorMapEntry color="#927D9E" quantity="39" label="Haplic Cryosols" opacity="0.7"/>
              <ColorMapEntry color="#FC9C4D" quantity="44" label="Haplic Ferralsols" opacity="0.7"/>
              <ColorMapEntry color="#E36940" quantity="45" label="Haplic Ferralsols (Rhodic)" opacity="0.7"/>
              <ColorMapEntry color="#FC934D" quantity="46" label="Haplic Ferralsols (Xanthic)" opacity="0.7"/>
              <ColorMapEntry color="#10A9E9" quantity="48" label="Haplic Fluvisols" opacity="0.7"/>
              <ColorMapEntry color="#0FC9E9" quantity="49" label="Haplic Fluvisols (Arenic)" opacity="0.7"/>
              <ColorMapEntry color="#0F8BE9" quantity="50" label="Haplic Fluvisols (Calcaric)" opacity="0.7"/>
              <ColorMapEntry color="#0F49E9" quantity="51" label="Haplic Fluvisols (Dystric)" opacity="0.7"/>
              <ColorMapEntry color="#0F8EE9" quantity="52" label="Haplic Fluvisols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#A087BF" quantity="54" label="Haplic Gleysols" opacity="0.7"/>
              <ColorMapEntry color="#793FBF" quantity="55" label="Haplic Gleysols (Dystric)" opacity="0.7"/>
              <ColorMapEntry color="#9E83BF" quantity="56" label="Haplic Gleysols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#FFF587" quantity="60" label="Haplic Gypsisols" opacity="0.7"/>
              <ColorMapEntry color="#D28769" quantity="67" label="Haplic Kastanozems" opacity="0.7"/>
              <ColorMapEntry color="#988F98" quantity="68" label="Haplic Leptosols" opacity="0.7"/>
              <ColorMapEntry color="#C9C9C9" quantity="69" label="Haplic Leptosols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#F7C7CA" quantity="73" label="Haplic Lixisols" opacity="0.7"/>
              <ColorMapEntry color="#F7B5E9" quantity="74" label="Haplic Lixisols (Chromic)" opacity="0.7"/>
              <ColorMapEntry color="#F795B6" quantity="75" label="Haplic Lixisols (Ferric)" opacity="0.7"/>
              <ColorMapEntry color="#F99491" quantity="79" label="Haplic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#F9919B" quantity="80" label="Haplic Luvisols (Chromic)" opacity="0.7"/>
              <ColorMapEntry color="#F99C92" quantity="81" label="Haplic Luvisols (Ferric)" opacity="0.7"/>
              <ColorMapEntry color="#CC7B68" quantity="86" label="Haplic Nitisols (Rhodic)" opacity="0.7"/>
              <ColorMapEntry color="#D5EADB" quantity="87" label="Haplic Phaeozems" opacity="0.7"/>
              <ColorMapEntry color="#F0966B" quantity="91" label="Haplic Planosols (Dystric)" opacity="0.7"/>
              <ColorMapEntry color="#F69C69" quantity="92" label="Haplic Planosols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#6DAA61" quantity="98" label="Haplic Podzols" opacity="0.7"/>
              <ColorMapEntry color="#FDD897" quantity="101" label="Haplic Regosols (Dystric)" opacity="0.7"/>
              <ColorMapEntry color="#FDDA9C" quantity="102" label="Haplic Regosols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#FFD49C" quantity="103" label="Haplic Regosols (Sodic)" opacity="0.7"/>
              <ColorMapEntry color="#E05BB4" quantity="106" label="Haplic Solonchaks" opacity="0.7"/>
              <ColorMapEntry color="#E05B7F" quantity="107" label="Haplic Solonchaks (Sodic)" opacity="0.7"/>
              <ColorMapEntry color="#F8AADE" quantity="110" label="Haplic Solonetz" opacity="0.7"/>
              <ColorMapEntry color="#8C7E79" quantity="113" label="Haplic Umbrisols" opacity="0.7"/>
              <ColorMapEntry color="#A88E99" quantity="116" label="Haplic Vertisols" opacity="0.7"/>
              <ColorMapEntry color="#A85C7D" quantity="117" label="Haplic Vertisols (Eutric)" opacity="0.7"/>
              <ColorMapEntry color="#3C5A6A" quantity="64" label="Hemic Histosols" opacity="0.7"/>
              <ColorMapEntry color="#FFD372" quantity="8" label="Histic Albeluvisols" opacity="0.7"/>
              <ColorMapEntry color="#FEEEC0" quantity="19" label="Hypoluvic Arenosols" opacity="0.7"/>
              <ColorMapEntry color="#FDE170" quantity="34" label="Leptic Cambisols" opacity="0.7"/>
              <ColorMapEntry color="#F99184" quantity="82" label="Leptic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#B4EADB" quantity="88" label="Leptic Phaeozems" opacity="0.7"/>
              <ColorMapEntry color="#FDD58E" quantity="104" label="Leptic Regosols" opacity="0.7"/>
              <ColorMapEntry color="#798C83" quantity="114" label="Leptic Umbrisols" opacity="0.7"/>
              <ColorMapEntry color="#B4BDB6" quantity="70" label="Lithic Leptosols" opacity="0.7"/>
              <ColorMapEntry color="#AB7269" quantity="96" label="Lixic Plinthosols" opacity="0.7"/>
              <ColorMapEntry color="#FFFC2B" quantity="23" label="Luvic Calcisols" opacity="0.7"/>
              <ColorMapEntry color="#E5D85C" quantity="38" label="Luvic Chernozems" opacity="0.7"/>
              <ColorMapEntry color="#C2EADB" quantity="89" label="Luvic Phaeozems" opacity="0.7"/>
              <ColorMapEntry color="#D29964" quantity="93" label="Luvic Planosols" opacity="0.7"/>
              <ColorMapEntry color="#73C3F4" quantity="112" label="Luvic Stagnosols" opacity="0.7"/>
              <ColorMapEntry color="#8759BF" quantity="57" label="Mollic Gleysols" opacity="0.7"/>
              <ColorMapEntry color="#9098A3" quantity="71" label="Mollic Leptosols" opacity="0.7"/>
              <ColorMapEntry color="#F8DFE6" quantity="111" label="Mollic Solonetz" opacity="0.7"/>
              <ColorMapEntry color="#A851A2" quantity="118" label="Mollic Vertisols" opacity="0.7"/>
              <ColorMapEntry color="#FFEE2B" quantity="24" label="Petric Calcisols" opacity="0.7"/>
              <ColorMapEntry color="#FAEBC6" quantity="42" label="Petric Durisols" opacity="0.7"/>
              <ColorMapEntry color="#EC6801" quantity="5" label="Plinthic Acrisols" opacity="0.7"/>
              <ColorMapEntry color="#FEF3C0" quantity="20" label="Protic Arenosols" opacity="0.7"/>
              <ColorMapEntry color="#A2A2A2" quantity="72" label="Rendzic Leptosols" opacity="0.7"/>
              <ColorMapEntry color="#7F807A" quantity="65" label="Sapric Histosols" opacity="0.7"/>
              <ColorMapEntry color="#D8955D" quantity="94" label="Solodic Planosols" opacity="0.7"/>
              <ColorMapEntry color="#F3B992" quantity="83" label="Stagnic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#7A7C88" quantity="40" label="Turbic Cryosols" opacity="0.7"/>
              <ColorMapEntry color="#FFD170" quantity="9" label="Umbric Albeluvisols" opacity="0.7"/>
              <ColorMapEntry color="#F59641" quantity="47" label="Umbric Ferralsols" opacity="0.7"/>
              <ColorMapEntry color="#936FBF" quantity="58" label="Umbric Gleysols" opacity="0.7"/>
              <ColorMapEntry color="#FDDF70" quantity="35" label="Vertic Cambisols" opacity="0.7"/>
              <ColorMapEntry color="#F99192" quantity="84" label="Vertic Luvisols" opacity="0.7"/>
              <ColorMapEntry color="#EC6806" quantity="6" label="Vetic Acrisols" opacity="0.7"/>
              <ColorMapEntry color="#FC5546" quantity="14" label="Vitric Andosols" opacity="0.7"/>
              <ColorMapEntry color="#9E8C7D" quantity="41" label="Vitric Cryosols" opacity="0.7"/>
            </ColorMap>
          </RasterSymbolizer>
        </Rule>
      </FeatureTypeStyle>
    </UserStyle>
  </NamedLayer>
</StyledLayerDescriptor>