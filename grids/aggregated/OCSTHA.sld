<?xml version="1.0" ?>
<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>OCSTHA</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                     <Geometry>
                      <ogc:PropertyName>Soil organic carbon stock in tonnes per ha</ogc:PropertyName>
                    </Geometry>
                        <ColorMap>
                            <ColorMapEntry color="#FFFFFF" quantity="-32768" label="NODATA" opacity="0.0"/>
                            <ColorMapEntry color="#0000ff" label="0" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#0050af" label="56" opacity="1.0" quantity="25"/>
                            <ColorMapEntry color="#00a15e" label="56" opacity="1.0" quantity="50"/>
                            <ColorMapEntry color="#00f20d" label="113" opacity="1.0" quantity="85"/>
                            <ColorMapEntry color="#43ff00" label="169" opacity="1.0" quantity="140"/>
                            <ColorMapEntry color="#94ff00" label="225" opacity="1.0" quantity="180"/>
                            <ColorMapEntry color="#e5ff00" label="281" opacity="1.0" quantity="220"/>
                            <ColorMapEntry color="#ffca00" label="338" opacity="1.0" quantity="280"/>
                            <ColorMapEntry color="#ff7900" label="394" opacity="1.0" quantity="320"/>
                            <ColorMapEntry color="#ff7900" label="394" opacity="1.0" quantity="400"/>
                            <ColorMapEntry color="#ff0000" label="850" opacity="1.0" quantity="850"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
