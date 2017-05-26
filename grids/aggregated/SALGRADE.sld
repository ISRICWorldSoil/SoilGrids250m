<?xml version="1.0" ?>
<StyledLayerDescriptor xmlns="http://www.opengis.net/sld" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml" xsi:schemaLocation="http://www.opengis.net/sld StyledLayerDescriptor.xsd" version="1.0.0">
    <UserLayer>
        <LayerFeatureConstraints>
            <FeatureTypeConstraint/>
        </LayerFeatureConstraints>
        <UserStyle>
            <Name>SLGWRB</Name>
            <Title/>
            <FeatureTypeStyle>
                <Name/>
                <Rule>
                    <RasterSymbolizer>
                     <Geometry>
                      <ogc:PropertyName>Saline / sodic soil grade based on WRB soil types and soil pH (SoilGrids250m)</ogc:PropertyName>
                    </Geometry>
                        <ColorMap>
                            <ColorMapEntry color="#f7f4f9" label="No saline/sodic properties" opacity="1.0" quantity="0"/>
                            <ColorMapEntry color="#d4b9da" label="Some sodic properties" opacity="1.0" quantity="1"/>
                            <ColorMapEntry color="#df65b0" label="Moderate sodic and solodic properties" opacity="1.0" quantity="2"/>
                            <ColorMapEntry color="#ce1256" label="Solonetz or Calcisols, sodic properties and/or ph &gt;8.1" opacity="1.0" quantity="3"/>
                            <ColorMapEntry color="#67001f" label="Solonetz, Solonchak, gypsic properties and/or ph &gt;8.5" opacity="1.0" quantity="4"/>
                            <ColorMapEntry color="#FFFFFF" quantity="255" label="NODATA" opacity="0.0"/>
                        </ColorMap>
                    </RasterSymbolizer>
                </Rule>
            </FeatureTypeStyle>
        </UserStyle>
    </UserLayer>
</StyledLayerDescriptor>
