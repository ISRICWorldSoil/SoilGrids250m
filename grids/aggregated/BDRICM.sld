<?xml version="1.0" ?>
<StyledLayerDescriptor version="1.0.0" xmlns="http://www.opengis.net/sld" xmlns:gml="http://www.opengis.net/gml" xmlns:ogc="http://www.opengis.net/ogc" xmlns:sld="http://www.opengis.net/sld">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>BDRICM</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                    <Geometry>
              <PropertyName>Depth to bedrock (R horizon) up to 200 cm</PropertyName>
              <Opacity>1</Opacity>
            </Geometry>
                        <ColorMap>
                        
                            <ColorMapEntry color="#b30000" label="0.0" opacity="0.7" quantity="0"/>
                            <ColorMapEntry color="#ce2a1d" label="29" opacity="0.7" quantity="28.5714"/>
                            <ColorMapEntry color="#e65338" label="57" opacity="0.7" quantity="57.1429"/>
                            <ColorMapEntry color="#f4794e" label="86" opacity="0.7" quantity="85.7143"/>
                            <ColorMapEntry color="#fc9f67" label="114" opacity="0.7" quantity="114.286"/>
                            <ColorMapEntry color="#fcc383" label="143" opacity="0.7" quantity="142.857"/>
                            <ColorMapEntry color="#fddbab" label="171" opacity="0.7" quantity="171.429"/>
                            <ColorMapEntry color="#fef0d9" label="200" opacity="0.7" quantity="200"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
